"""Meaningful exhaustive and seeded checks against a physical-state oracle."""

import ast
import json
from pathlib import Path
import random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from boolean_certificates import (BooleanModel, LocalRule, all_minimum_copy_sets,
                                  cube_row_table, load_json_model,
                                  local_rule_from_callable, physical_residual_model,
                                  pure_copy_certificate, support_map, support_trace)
from independent_oracle import all_subsets, classify
from extract_models import safe_arithmetic


def bnet_evaluate(formula, state):
    """Independent parser for emitted BNET formulas, using a restricted AST."""
    expression = formula.replace("!", " not ").replace("&", " and ").replace("|", " or ").strip()
    tree = ast.parse(expression, mode="eval")

    def evaluate(node):
        if isinstance(node, ast.Expression):
            return evaluate(node.body)
        if isinstance(node, ast.Name) and node.id.startswith("n") and node.id[1:].isdigit():
            return bool(state[int(node.id[1:]) - 1])
        if isinstance(node, ast.Constant) and node.value in (0, 1):
            return bool(node.value)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
            return not evaluate(node.operand)
        if isinstance(node, ast.BoolOp) and isinstance(node.op, ast.And):
            return all(evaluate(x) for x in node.values)
        if isinstance(node, ast.BoolOp) and isinstance(node.op, ast.Or):
            return any(evaluate(x) for x in node.values)
        raise ValueError("Unsupported BNET syntax")
    return int(evaluate(tree))


def validate_model(model, counts):
    before = repr(model)
    universe = frozenset(range(model.n))
    zero = model.relabeled()
    # Every support cube and output row is compared with direct global-state
    # evaluation, including correlations from repeated physical references.
    for support in all_subsets(universe):
        direct = set()
        for active in all_subsets(support):
            physical = tuple(model.target[i] ^ int(i in active) for i in range(model.n))
            output = model.step(physical)
            direct.update(i for i in universe if output[i] != model.target[i])
        if set(support_map(zero, support)) != direct:
            raise AssertionError("Cube support disagrees with independent physical enumeration")
        counts["support_cubes"] += 1
    for baseline in all_subsets(universe):
        trace = support_trace(model, baseline)
        residual = trace[-1]
        residual_model, node_order = physical_residual_model(model, residual)
        location = {node: index for index, node in enumerate(node_order)}
        optimum = None
        successful_sets = []
        for additional in all_subsets(universe - baseline):
            original = classify(model.rules, model.target, baseline | additional)
            if residual_model is None:
                reduced_success = True
            else:
                reduced_pins = [location[p] for p in additional & residual]
                reduced_success = classify(residual_model.rules, residual_model.target,
                                           reduced_pins)["globally_convergent"]
            if original["globally_convergent"] != reduced_success:
                raise AssertionError("Physical residual success equivalence failed")
            if original["globally_convergent"]:
                successful_sets.append(additional)
                optimum = len(additional) if optimum is None else min(optimum, len(additional))
            counts["baseline_additional_pairs"] += 1
        certificate = pure_copy_certificate(model, baseline)
        if certificate is not None:
            if certificate["conditional_optimum"] != optimum:
                raise AssertionError("Pure-copy cost differs from full-state exact optimum")
            predicted = {frozenset(k) for k in all_minimum_copy_sets(certificate)}
            actual = {k for k in successful_sets if len(k) == optimum}
            if predicted != actual:
                raise AssertionError("Pure-copy minimum sets differ from exhaustive optimum sets")
            counts["pure_copy_certificates"] += 1
    counts["networks"] += 1
    if repr(model) != before:
        raise AssertionError("Certificate routines mutated an input model")


def validate_all(models_directory):
    counts = {"networks": 0, "support_cubes": 0, "baseline_additional_pairs": 0,
              "pure_copy_certificates": 0}
    # All two-node Boolean maps compatible with each of the four possible
    # targets, allowing both physical coordinates in each truth table.
    for target_code in range(4):
        target = tuple((target_code >> i) & 1 for i in range(2))
        admissible = []
        for i in range(2):
            admissible.append([tuple((table_code >> k) & 1 for k in range(4))
                               for table_code in range(16)
                               if ((table_code >> target_code) & 1) == target[i]])
        for first in admissible[0]:
            for second in admissible[1]:
                validate_model(BooleanModel((LocalRule((0, 1), first),
                                              LocalRule((0, 1), second)), target), counts)
    exhaustive_two_node_count = counts["networks"]
    rng = random.Random(2026101401)
    for n, repetitions in ((3, 32), (4, 32), (5, 8)):
        for _ in range(repetitions):
            target = tuple(rng.randrange(2) for _ in range(n))
            rows = []
            for i in range(n):
                arity = rng.randrange(n + 1)
                parents = tuple(rng.sample(range(n), arity))
                entries = [rng.randrange(2) for _ in range(1 << arity)]
                target_corner = sum(target[p] << j for j, p in enumerate(parents))
                entries[target_corner] = target[i]
                rows.append(LocalRule(parents, tuple(entries)))
            validate_model(BooleanModel(tuple(rows), target), counts)
    # Repeated references must collapse to one physical variable. These
    # expressions would produce spurious supports under independent-slot rules.
    repeated_zero = local_rule_from_callable((0, 0), lambda a, b: a and not b)
    repeated_copy = local_rule_from_callable((0, 0, 1, 1),
                                            lambda a, b, c, d: (a and b) or (c and not d))
    assert repeated_zero == LocalRule((0,), (0, 0))
    assert repeated_copy == LocalRule((0, 1), (0, 1, 0, 1))
    repeated_model = BooleanModel((repeated_zero, repeated_copy), (0, 0))
    validate_model(repeated_model, counts)
    # Global success is not monotone under additional pins in general.
    nonmonotone = BooleanModel((local_rule_from_callable((0, 1), lambda a, b: a and not b),
                                LocalRule((0,), (0, 1))), (0, 0))
    unpinned = classify(nonmonotone.rules, nonmonotone.target)
    pinned_second = classify(nonmonotone.rules, nonmonotone.target, {1})
    assert unpinned["globally_convergent"] and not pinned_second["globally_convergent"]
    # Input validation is part of the semantics contract.
    invalid = 0
    for thunk in (lambda: LocalRule((0, 0), (0, 0, 0, 1)),
                  lambda: LocalRule((0,), (0, 1, 0)),
                  lambda: LocalRule((0,), (0, 2)),
                  lambda: BooleanModel((LocalRule((0,), (1, 0)),), (0,)),
                  lambda: BooleanModel((LocalRule((1,), (0, 1)),), (0,))):
        try:
            thunk()
        except ValueError:
            invalid += 1
        else:
            raise AssertionError("Invalid input was accepted")
    # Exported BNET formulas are checked independently against every local corner.
    exported_checks = {}
    polynomial_checks = {}
    for file in sorted(Path(models_directory).glob("*.json")):
        document = json.loads(file.read_text())
        model = load_json_model(document)
        count = 0
        polynomial_count = 0
        bnet_rows = [line for line in file.with_suffix(".bnet").read_text().splitlines()
                     if line and not line.startswith("#") and not line.startswith("targets,")]
        if len(bnet_rows) != model.n:
            raise AssertionError("BNET row count mismatch")
        for i, (row, rule, line) in enumerate(zip(document["rules"], model.rules, bnet_rows)):
            label, formula = line.split(",", 1)
            assert label.strip() == f"n{i+1:02d}" and formula.strip() == row["bnet_formula"]
            for code, expected in enumerate(rule.table):
                state = [0] * model.n
                for j, p in enumerate(rule.parents):
                    state[p] = (code >> j) & 1
                if bnet_evaluate(formula, state) != expected:
                    raise AssertionError("BNET formula disagrees with explicit physical truth table")
                count += 1
            source_expression = row.get("source_polynomial_expression", row.get("source_expression"))
            if source_expression is not None:
                expression = ast.parse(source_expression, mode="eval").body
                source_parents = sorted({ast.literal_eval(node.slice) for node in ast.walk(expression)
                                         if isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name)
                                         and node.value.id == "S"})
                for code in range(1 << len(source_parents)):
                    values = {p: (code >> j) & 1 for j, p in enumerate(source_parents)}
                    state = [0] * model.n
                    for p, value in values.items():
                        state[p] = value
                    expected = safe_arithmetic(expression, values)
                    if expected != bnet_evaluate(formula, state):
                        raise AssertionError("Original polynomial corner differs from explicit BNET rule")
                    polynomial_count += 1
        exported_checks[document["model_id"]] = count
        polynomial_checks[document["model_id"]] = polynomial_count
    return {"status": "PASS", "scope": "Finite exhaustive validation and deterministic seeded implementation checks.",
            "seed": 2026101401, "exhaustive_two_node_target_map_pairs": exhaustive_two_node_count,
            "seeded_three_four_five_node_models": 72, "repeated_reference_models": 1,
            "counts": counts, "invalid_inputs_rejected": invalid,
            "input_models_unchanged": True,
            "nonmonotone_feasibility_counterexample_confirmed": True,
            "exported_bnet_local_corners_checked": exported_checks,
            "original_source_polynomial_corners_checked": polynomial_checks}


if __name__ == "__main__":
    print(json.dumps(validate_all(Path(__file__).resolve().parents[1] / "models"), indent=2))
