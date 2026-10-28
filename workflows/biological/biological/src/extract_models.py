"""Read local source files as AST data, never import/execute the source scripts.

This optional extraction utility recreates the supplied JSON/BNET model files.
The certificate reproducer uses the explicit model files and bundled sources.
"""

import argparse
import ast
import hashlib
import json
from pathlib import Path

from boolean_certificates import BooleanModel, LocalRule, local_rule_from_callable


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def template(kind, values):
    if kind == "copy":
        return values[0]
    if kind == "or_not":
        return int(values[0] or not values[1])
    if kind == "or_and_or":
        return int(values[0] or (values[1] and (values[2] or values[3])))
    if kind == "notandnot":
        return int(not values[0] and not values[1])
    if kind == "or":
        return int(values[0] or values[1])
    if kind == "not":
        return 1 - values[0]
    if kind == "notand":
        return int(not values[0] and values[1])
    raise ValueError("Unrecognized source rule template")


def template_formula(kind, slots):
    s = [f"n{i+1:02d}" for i in slots]
    return {
        "copy": lambda: s[0],
        "or_not": lambda: f"({s[0]} | !{s[1]})",
        "or_and_or": lambda: f"({s[0]} | ({s[1]} & ({s[2]} | {s[3]})))",
        "notandnot": lambda: f"(!{s[0]} & !{s[1]})",
        "or": lambda: f"({s[0]} | {s[1]})",
        "not": lambda: f"!{s[0]}",
        "notand": lambda: f"(!{s[0]} & {s[1]})",
    }[kind]()


def truth_table_dnf(rule):
    if all(rule.table):
        return "1"
    if not any(rule.table):
        return "0"
    terms = []
    for code, value in enumerate(rule.table):
        if value:
            terms.append("(" + " & ".join(("" if (code >> j) & 1 else "!")
                                          + f"n{p+1:02d}"
                                          for j, p in enumerate(rule.parents)) + ")")
    return " | ".join(terms)


def model_document(model, rows, provenance, key):
    return {
        "schema": "BooleanPhysicalTruthTables_v1",
        "model_id": key,
        "n": model.n,
        "semantics": "Synchronous deterministic Boolean updates; static permanent target-value pins from time zero; every compatible initial state; unit costs; all nodes eligible.",
        "index_convention": "JSON parents and engine use zero-based source indices. BNET labels n01,... are one-based. Truth-table index bit j is the value of parents_zero_based[j].",
        "provenance": provenance,
        "target": list(model.target),
        "target_node_order": [f"n{i+1:02d}" for i in range(model.n)],
        "target_integer_lsb_node1": sum(v << i for i, v in enumerate(model.target)),
        "rules": rows,
    }


def extract29(template_source, polynomial_source=None):
    content = Path(template_source).read_bytes()
    source = content.decode("utf-8")
    tree = ast.parse(source)
    rule_node = next(n.value for n in tree.body if isinstance(n, ast.AnnAssign)
                     and isinstance(n.target, ast.Name) and n.target.id == "TRUE_RULES")
    target_node = next(n.value for n in tree.body if isinstance(n, ast.Assign)
                       and any(isinstance(t, ast.Name) and t.id == "TARGET" for t in n.targets))
    if not (isinstance(rule_node, ast.List) and isinstance(target_node, ast.Call)):
        raise ValueError("Unexpected 29-node AST layout")
    target = tuple(int(x) for x in ast.literal_eval(target_node.args[0]))
    rules, rows = [], []
    for index, call in enumerate(rule_node.elts):
        if not isinstance(call, ast.Call) or not isinstance(call.func, ast.Name) or call.func.id != "Rule":
            raise ValueError("Unexpected 29-node rule constructor")
        kind, slots = ast.literal_eval(call.args[0]), ast.literal_eval(call.args[1])
        rule = local_rule_from_callable(slots, lambda *values, kind=kind: template(kind, values))
        rules.append(rule)
        rows.append({"node_zero_based": index, "node_one_based": index + 1,
                     "node_label": f"n{index+1:02d}", "source_template": kind,
                     "source_references_zero_based": slots,
                     "source_line": call.lineno,
                     "parents_zero_based": list(rule.parents),
                     "physical_truth_table": list(rule.table),
                     "bnet_formula": template_formula(kind, slots)})
    model = BooleanModel(tuple(rules), target)
    if model.n != 29:
        raise ValueError("Unexpected network size")
    document = model_document(model, rows, {
        "source_path": "provenance/source_files/" + Path(template_source).name,
        "source_file_sha256": sha256(content),
        "authority": "Rule templates explicitly described in the supplied source as inferred from the released 29-node QL script. This artifact certifies that coded Boolean implementation, not fidelity to a published biological reconstruction.",
        "source_extraction": "AST literal extraction only; the source script was not imported or executed.",
        "biological_gene_labels": "Not asserted; nodes retain neutral source-order indices.",
    }, "tlgl_coded_29")
    if polynomial_source is not None:
        content = Path(polynomial_source).read_bytes()
        source = content.decode("utf-8-sig")
        tree = ast.parse(source)
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "calculate_H")
        assignments = {}
        for statement in ast.walk(function):
            if not isinstance(statement, ast.Assign) or len(statement.targets) != 1:
                continue
            dest = statement.targets[0]
            if (isinstance(dest, ast.Subscript) and isinstance(dest.value, ast.Name)
                    and dest.value.id == "new_S" and isinstance(dest.slice, ast.Constant)
                    and type(dest.slice.value) is int):
                if dest.slice.value in assignments:
                    raise ValueError("Repeated literal polynomial row assignment")
                assignments[dest.slice.value] = statement
        if set(assignments) != set(range(29)):
            raise ValueError("Incomplete 29-node polynomial assignment block")
        corner_count = 0
        for index in range(29):
            expression = assignments[index].value
            parents = sorted({ast.literal_eval(n.slice) for n in ast.walk(expression)
                              if isinstance(n, ast.Subscript) and isinstance(n.value, ast.Name)
                              and n.value.id == "S"})
            row = document["rules"][index]
            for code in range(1 << len(parents)):
                values = {p: (code >> j) & 1 for j, p in enumerate(parents)}
                physical = [0] * 29
                for p, value in values.items():
                    physical[p] = value
                if safe_arithmetic(expression, values) != model.rules[index].evaluate(physical):
                    raise AssertionError("Inferred template differs from the original coded polynomial")
                corner_count += 1
            row["source_polynomial_expression"] = ast.get_source_segment(source, expression)
            row["source_polynomial_references_zero_based"] = parents
        target_node = next(n.value for n in tree.body if isinstance(n, ast.Assign)
                           and any(isinstance(t, ast.Name) and t.id == "attractor" for t in n.targets))
        if tuple(int(x) for x in ast.literal_eval(target_node)) != model.target:
            raise AssertionError("The original coded target differs from the template target")
        document["provenance"]["original_polynomial_crosscheck"] = {
            "source_file_name": Path(polynomial_source).name,
            "source_file_sha256": sha256(content),
            "source_path": "provenance/source_files/29nodeQL.py",
            "all_source_local_corners_checked": corner_count,
            "status": "PASS",
            "corner_count_note": "The original polynomial at one-based node 18 lists three source variables but simplifies to NOT n17. Its eight source corners reduce to two physical truth-table corners, explaining 88 source corners versus 82 exported template corners.",
        }
    return document


def safe_arithmetic(node, values):
    """Integer-only interpreter for the source polynomial assignment ASTs."""
    if isinstance(node, ast.Constant) and type(node.value) is int:
        return node.value
    if isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name) and node.value.id == "S":
        index = ast.literal_eval(node.slice)
        if type(index) is not int:
            raise ValueError("Nonliteral source index")
        return values[index]
    if isinstance(node, ast.BinOp):
        left, right = safe_arithmetic(node.left, values), safe_arithmetic(node.right, values)
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -safe_arithmetic(node.operand, values)
    raise ValueError(f"Unsupported source syntax: {type(node).__name__}")


def extract17(source_path):
    content = Path(source_path).read_bytes()
    source = content.decode("utf-8-sig")
    tree = ast.parse(source)
    target_value = next(n.value for n in tree.body if isinstance(n, ast.Assign)
                        and any(isinstance(t, ast.Name) and t.id == "attractor" for t in n.targets))
    target = tuple(int(x) for x in ast.literal_eval(target_value))
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "calculate_H")
    assignments = {}
    for statement in ast.walk(function):
        if not isinstance(statement, ast.Assign) or len(statement.targets) != 1:
            continue
        dest = statement.targets[0]
        if (isinstance(dest, ast.Subscript) and isinstance(dest.value, ast.Name)
                and dest.value.id == "new_S" and isinstance(dest.slice, ast.Constant)
                and type(dest.slice.value) is int):
            if dest.slice.value in assignments:
                raise ValueError("Repeated literal row assignment")
            assignments[dest.slice.value] = statement
    if set(assignments) != set(range(17)):
        raise ValueError("Incomplete 17-node assignment block")
    rules, rows = [], []
    for index in range(17):
        statement = assignments[index]
        parents = tuple(sorted({ast.literal_eval(n.slice) for n in ast.walk(statement.value)
                                if isinstance(n, ast.Subscript) and isinstance(n.value, ast.Name)
                                and n.value.id == "S"}))
        table = tuple(safe_arithmetic(statement.value, {p: (code >> j) & 1
                                                       for j, p in enumerate(parents)})
                      for code in range(1 << len(parents)))
        rule = LocalRule(parents, table)
        rules.append(rule)
        rows.append({"node_zero_based": index, "node_one_based": index + 1,
                     "node_label": f"n{index+1:02d}",
                     "source_expression": ast.get_source_segment(source, statement.value),
                     "source_lines": [statement.lineno, statement.end_lineno],
                     "parents_zero_based": list(parents),
                     "physical_truth_table": list(table),
                     "bnet_formula": truth_table_dnf(rule)})
    return model_document(BooleanModel(tuple(rules), target), rows, {
        "source_file_name": Path(source_path).name,
        "source_file_sha256": sha256(content),
        "source_path": "provenance/source_files/17nodeQL.py",
        "authority": "Boolean corners of the 17-node coded polynomial assignments before the automatic deterministic-coordinate freezing loop; target copied from its attractor literal. No biological rule correction is inferred.",
        "source_extraction": "AST parsing and an integer arithmetic whitelist; the source script was not imported or executed.",
        "update_semantics_change": "The original learning/IBMF automatic-freezing procedure is not run. These certificates concern the explicit synchronous Boolean map with only the specified permanent pins.",
        "biological_gene_labels": "Not asserted; nodes retain neutral source-order indices.",
    }, "drosophila_coded_17")


def write_model(document, directory):
    path = Path(directory)
    stem = document["model_id"]
    json_path, bnet_path = path / (stem + ".json"), path / (stem + ".bnet")
    if json_path.exists() or bnet_path.exists():
        raise FileExistsError("Refusing to overwrite an existing extracted model")
    json_path.write_text(json.dumps(document, indent=2) + "\n")
    lines = ["# Exact physical Boolean rules; neutral one-based source labels.",
             "# Permanent target pins and target values are specified in the matching JSON.",
             "targets, factors"]
    lines.extend(f'{r["node_label"]}, {r["bnet_formula"]}' for r in document["rules"])
    bnet_path.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source29-template", required=True, help="Bundled T-LGL rule-template source")
    parser.add_argument("--source17", required=True)
    parser.add_argument("--source29", help="Optional frozen original 29nodeQL.py for independent polynomial-corner comparison")
    parser.add_argument("--output", required=True, help="Existing empty output directory")
    args = parser.parse_args()
    output = Path(args.output)
    if not output.is_dir() or any(output.iterdir()):
        parser.error("--output must be an existing empty directory")
    write_model(extract29(args.source29_template, args.source29), output)
    write_model(extract17(args.source17), output)


if __name__ == "__main__":
    main()
