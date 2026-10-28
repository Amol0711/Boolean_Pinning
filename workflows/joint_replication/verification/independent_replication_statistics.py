#!/usr/bin/env python3
"""Independent all-draw statistics for the prespecified J1 replication J1 replication.

No production modules are imported. All moment numerators are exact integer
sums; law parameters are rederived with fractions. Diagnostics against the
declared law never estimate or replace the reference Poisson parameter.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist

from scipy.special import betaincinv
from scipy.stats import chi2


Z95 = NormalDist().inv_cdf(.975)
RULES = ((0, "ZERO"), (5, "COPY1"), (6, "COPY2"), (4, "AND"),
         (7, "OR"), (3, "XOR"))
DECLARED_WEIGHTS = {1: [30, 0, 0, 0, 15, 45, 10, 0],
                    3: [60, 0, 0, 0, 10, 10, 20, 0]}


class Audit:
    def __init__(self):
        self.checks = Counter()
        self.failures = []
        self.hashes = {}

    def check(self, condition, category, label):
        self.checks[category] += 1
        if not condition:
            self.failures.append(label)

    def hash(self, path, root):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        try:
            label = str(path.relative_to(root))
        except ValueError:
            label = str(path)
        self.hashes[label] = digest
        return digest


def read_rows(path):
    with gzip.open(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream))


def cp_interval(successes, trials, alpha=.05):
    """Equal-tailed Clopper--Pearson interval; no asymptotic approximation."""
    lower = 0. if successes == 0 else float(betaincinv(successes, trials-successes+1, alpha/2))
    upper = 1. if successes == trials else float(betaincinv(successes+1, trials-successes, 1-alpha/2))
    return [lower, upper]


def endpoint_moment(values, power):
    n = len(values)
    total = sum(value**power for value in values)
    square_total = sum(value**(2*power) for value in values)
    mean = total/n
    se = math.sqrt((n*square_total-total*total)/(n*n*(n-1)))
    return {"sum": total, "sum_of_squared_variables": square_total,
            "average": mean, "se": se}


def summarize(rows, label):
    n = len(rows)
    lower = [int(row["cost_lower"]) for row in rows]
    upper = [int(row["cost_upper"]) for row in rows]
    out = {"label": label, "family": "J1", "n": 25600, "draws": n,
           "unresolved": sum(a != b for a, b in zip(lower, upper)),
           "histogram_lower": dict(sorted(Counter(lower).items())),
           "histogram_upper": dict(sorted(Counter(upper).items()))}
    for name, power in (("mean", 1), ("m2", 2)):
        lo, hi = endpoint_moment(lower, power), endpoint_moment(upper, power)
        out[name] = {"lower": lo["average"], "upper": hi["average"],
                     "se_lower": lo["se"], "se_upper": hi["se"],
                     "normal95_envelope": [max(0., lo["average"]-Z95*lo["se"]),
                                            hi["average"]+Z95*hi["se"]],
                     "integer_sums_lower": {k: v for k, v in lo.items() if k.startswith("sum")},
                     "integer_sums_upper": {k: v for k, v in hi.items() if k.startswith("sum")}}
    definite, possible = upper.count(0), lower.count(0)
    out["P0"] = {"definite_zero_count": definite, "possible_zero_count": possible,
                 "lower": definite/n, "upper": possible/n,
                 "cp95_envelope": [cp_interval(definite, n)[0], cp_interval(possible, n)[1]]}
    for flag in ("G", "G_prime"):
        failures = sum(int(row[flag]) == 0 for row in rows)
        out[flag+"_failures"] = failures
        out[flag+"_failures_cp95"] = cp_interval(failures, n)
    stops = Counter()
    for row in rows:
        stops.update(json.loads(row["stop_summary"]))
    out.update({"max_residual": max(int(row["residual_size"]) for row in rows),
                "max_complex": max(int(row["largest_complex"]) for row in rows),
                "total_solver_work": sum(int(row["solver_work"]) for row in rows),
                "stop_reason_counts": dict(sorted(stops.items())),
                "exact_reason_counts": dict(sorted(Counter(row["exact_reason"] for row in rows).items()))})
    return out


def derive_law(protocol, original, audit):
    family = protocol["families"]["J1"]
    base_family = original["families"]["J1"]
    lam = rho = mass = degree_average = Fraction(0)
    for index, typ in enumerate(family["types"]):
        degree = typ["degree"]
        fraction = Fraction(*typ["fraction"])
        weights = typ["weights"]
        denominator = typ["probability_denominator"]
        audit.check(weights == DECLARED_WEIGHTS.get(degree) and denominator == 100 and fraction == Fraction(1, 2),
                    "declared_J1_rule_law", str(index))
        audit.check(typ == base_family["types"][index], "unchanged_original_rule_law", str(index))
        p = [Fraction(w, denominator) for w in weights]
        q1 = sum((v for code, v in enumerate(p) if code & 1), Fraction(0))
        q2 = sum((v for code, v in enumerate(p) if code & 2), Fraction(0))
        lam += fraction*degree*(q1+q2)/2
        rho += fraction*degree*(1-p[0])/2
        mass += fraction
        degree_average += fraction*degree
    audit.check(mass == 1 and degree_average == 2 and lam == Fraction(29, 80) and rho == Fraction(19, 40),
                "exact_rational_parameters", "J1 rational derivation")
    audit.check(family["parameters_rational"] == base_family["parameters_rational"],
                "unchanged_reference", "protocol rational parameters")
    expected_parameters = {"lambda": [29, 80], "rho": [19, 40], "contraction_coefficient": [19, 40]}
    audit.check(family["parameters_rational"] == expected_parameters, "exact_rational_parameters", "protocol rational declarations")
    mu = -math.log1p(-float(lam))
    audit.check(abs(family["reference_mu"]-mu) < 2e-15 and family["reference_mu"] == base_family["reference_mu"],
                "unchanged_reference", "fixed mu")
    return {"lambda_rational": [29, 80], "rho_rational": [19, 40],
            "lambda": float(lam), "rho": float(rho), "mu": mu,
            "reference_m2": mu+mu*mu, "reference_p0": float(1-lam),
            "refitted": False}


def check_rows(rows, protocol, audit, baseline=False):
    expected = 2000 if baseline else 20000
    prefix = "original" if baseline else "replication"
    audit.check(len(rows) == expected and [int(row["replicate"]) for row in rows] == list(range(expected)),
                "complete_unique_ordered_draws", prefix)
    for row in rows:
        label = prefix+"/"+row["replicate"]
        lo, hi, residual = (int(row[k]) for k in ("cost_lower", "cost_upper", "residual_size"))
        audit.check(row["family"] == "J1" and int(row["n"]) == 25600,
                    "row_identity", label)
        audit.check(0 <= lo <= hi <= residual <= 25600, "cost_interval", label)
        g, gp = int(row["G"]), int(row["G_prime"])
        audit.check(g in (0, 1) and gp in (0, 1) and g <= gp and g == (int(row["double_rows"]) == 0),
                    "certificate_accounting", label)
        audit.check((gp and lo == hi == int(row["certificate_cycle_count"])) or
                    (not gp and row["certificate_cycle_count"] == ""), "certificate_cost", label)
        pins = [int(v) for v in row["pins"].split(";") if v]
        audit.check(len(pins) == len(set(pins)) == hi and all(0 <= v < 25600 for v in pins),
                    "witness_cardinality", label)
        reason = "original_G" if g else "weaker_G_prime" if gp else "SCC_exact" if lo == hi else "bounded_unresolved"
        audit.check(row["exact_reason"] == reason, "exact_reason", label)
        if not baseline:
            expected_entropy = [protocol["master_seed"], protocol["families"]["J1"]["id"], 25600, int(row["replicate"])]
            audit.check(json.loads(row["seed_entropy"]) == expected_entropy, "seed_entropy", label)
            for degree in (1, 3):
                counts = json.loads(row[f"rule_counts_d{degree}"])
                audit.check(len(counts) == 8 and all(type(v) is int and v >= 0 for v in counts)
                            and sum(counts) == 12800, "rule_count_accounting", label+f"/d{degree}")
                audit.check(all(counts[c] == 0 for c, w in enumerate(DECLARED_WEIGHTS[degree]) if w == 0),
                            "impossible_rule_categories", label+f"/d{degree}")


def detail_accounting(path, rows, audit):
    seen = 0
    with gzip.open(path, "rt") as stream:
        for seen, line in enumerate(stream, 1):
            detail = json.loads(line)
            if seen > len(rows):
                audit.check(False, "detail_accounting", "extra detail row")
                continue
            row = rows[seen-1]
            label = row["replicate"]
            audit.check(detail["family"] == row["family"] and detail["n"] == int(row["n"])
                        and detail["replicate"] == int(row["replicate"])
                        and detail["network_sha256"] == row["network_sha256"], "detail_identity", label)
            audit.check(detail["lower"] == int(row["cost_lower"]) and detail["upper"] == int(row["cost_upper"])
                        and detail["G"] == bool(int(row["G"])) and detail["G_prime"] == bool(int(row["G_prime"]))
                        and detail["exact_reason"] == row["exact_reason"]
                        and detail["pins"] == [int(v) for v in row["pins"].split(";") if v]
                        and len(detail["residual_rows"]) == int(row["residual_size"]),
                        "detail_cost_accounting", label)
            for key in ("seed_entropy", "rule_counts_d1", "rule_counts_d3", "component_sizes", "cyclic_component_sizes"):
                audit.check(detail[key] == json.loads(row[key]), "detail_extended_accounting", label+"/"+key)
    audit.check(seen == len(rows), "detail_accounting", "complete detail row count")


def rule_diagnostics(rows):
    summaries = []
    per_degree = []
    for degree in (1, 3):
        counts = [sum(json.loads(row[f"rule_counts_d{degree}"])[code] for row in rows) for code in range(8)]
        total = sum(counts)
        statistic = sum((counts[code]-total*w/100)**2/(total*w/100)
                        for code, w in enumerate(DECLARED_WEIGHTS[degree]) if w)
        per_degree.append({"degree": degree, "trials": total, "all_corner_counts": counts,
                           "pearson_chi_square": statistic, "degrees_of_freedom": 3,
                           "diagnostic_p_value": float(chi2.sf(statistic, 3))})
        for code, name in RULES:
            declared = DECLARED_WEIGHTS[degree][code]/100
            count = counts[code]
            pointwise = cp_interval(count, total) if declared else [0., 0.]
            simultaneous = cp_interval(count, total, .05/8) if declared else [0., 0.]
            summaries.append({"degree": degree, "rule": name, "corner_code": code,
                              "count": count, "total": total, "declared_probability": declared,
                              "observed_probability": count/total,
                              "expected_count": total*declared,
                              "standardized_count_residual": (count-total*declared)/math.sqrt(total*declared*(1-declared)) if declared else None,
                              "pointwise_cp95_lower": pointwise[0], "pointwise_cp95_upper": pointwise[1],
                              "simultaneous_cp95_lower": simultaneous[0], "simultaneous_cp95_upper": simultaneous[1],
                              "impossible_by_design": declared == 0,
                              "simultaneous_interval_contains_design_probability": simultaneous[0] <= declared <= simultaneous[1]})
    joint = sum(item["pearson_chi_square"] for item in per_degree)
    return {"rows": summaries, "per_degree": per_degree,
            "joint_pearson_chi_square": joint, "joint_degrees_of_freedom": 6,
            "joint_diagnostic_p_value": float(chi2.sf(joint, 6)),
            "scope": "Two independent multinomial samples with fixed class totals and fixed probabilities. Pearson p-values are asymptotic diagnostics, not correctness gates. Eight positive-probability categories share Bonferroni 95% simultaneous Clopper-Pearson coverage. Impossible categories are structural zeros, not estimated probabilities; all eight corner codes are counted."}


def differences(original, replication):
    result = {"direction": "replication minus original", "method": "Independent prespecified RNG streams; approximate pointwise normal uncertainty for finite-size population differences. The contrast is exploratory because J1 was selected after examining the original sample. No equality to the asymptotic benchmark is assumed."}
    for name in ("mean", "m2"):
        a, b = original[name], replication[name]
        lo, hi = b["lower"]-a["upper"], b["upper"]-a["lower"]
        se_lo = math.hypot(b["se_lower"], a["se_upper"])
        se_hi = math.hypot(b["se_upper"], a["se_lower"])
        result[name] = {"lower": lo, "upper": hi, "se_lower": se_lo, "se_upper": se_hi,
                        "normal95_envelope": [lo-Z95*se_lo, hi+Z95*se_hi],
                        "standardized_difference": lo/se_lo if lo == hi and se_lo else None}
    a, b = original["P0"], replication["P0"]
    low, high = b["lower"]-a["upper"], b["upper"]-a["lower"]
    se_lo = math.sqrt(b["lower"]*(1-b["lower"])/replication["draws"] + a["upper"]*(1-a["upper"])/original["draws"])
    se_hi = math.sqrt(b["upper"]*(1-b["upper"])/replication["draws"] + a["lower"]*(1-a["lower"])/original["draws"])
    result["P0"] = {"lower": low, "upper": high, "se_lower": se_lo, "se_upper": se_hi,
                    "normal95_envelope": [low-Z95*se_lo, high+Z95*se_hi]}
    return result


def compare_production(summary, expected, law, audit):
    targets = {k: expected[k] for k in ("family", "n", "draws", "unresolved", "G_failures", "G_prime_failures", "max_residual", "max_complex", "total_solver_work", "stop_reason_counts", "exact_reason_counts", "P0", "G_failures_cp95", "G_prime_failures_cp95")}
    targets.update({"lambda": law["lambda"], "rho": law["rho"], "contraction_coefficient": law["rho"], "mu": law["mu"], "reference_m2": law["reference_m2"], "reference_p0": law["reference_p0"]})
    for name in ("mean", "m2"):
        targets[name] = {k: expected[name][k] for k in ("lower", "upper", "se_lower", "se_upper", "normal95_envelope")}
    def compare(wanted, observed, label):
        if isinstance(wanted, dict):
            audit.check(isinstance(observed, dict), "production_summary_type", label)
            if isinstance(observed, dict):
                for key, value in wanted.items():
                    compare(value, observed.get(key), label+"/"+key)
        elif isinstance(wanted, list):
            audit.check(isinstance(observed, list) and len(observed) == len(wanted), "production_summary_type", label)
            if isinstance(observed, list):
                for index, (a, b) in enumerate(zip(wanted, observed)):
                    compare(a, b, label+f"/{index}")
        elif isinstance(wanted, float):
            audit.check(isinstance(observed, (int, float)) and math.isclose(wanted, observed, rel_tol=2e-10, abs_tol=2e-12), "production_summary_float", label)
        else:
            audit.check(wanted == observed, "production_summary_exact", label)
    compare(targets, summary, "summary")


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def compare_summary_csv(path, summary, law, audit):
    expected = {k: summary[k] for k in ("family", "n", "draws", "unresolved", "G_failures", "G_prime_failures", "max_residual", "max_complex", "total_solver_work")}
    expected.update({"lambda": law["lambda"], "rho": law["rho"], "contraction_coefficient": law["rho"],
                     "mu": law["mu"], "reference_m2": law["reference_m2"], "reference_p0": law["reference_p0"]})
    for label in ("mean", "m2"):
        for field in ("lower", "upper", "se_lower", "se_upper"):
            expected[label+"_"+field] = summary[label][field]
        expected[label+"_normal95_lo"], expected[label+"_normal95_hi"] = summary[label]["normal95_envelope"]
    for field in ("definite_zero_count", "possible_zero_count", "lower", "upper"):
        expected["P0_"+field] = summary["P0"][field]
    expected["P0_cp95_lo"], expected["P0_cp95_hi"] = summary["P0"]["cp95_envelope"]
    for label in ("G_failures", "G_prime_failures"):
        expected[label+"_cp95_lo"], expected[label+"_cp95_hi"] = summary[label+"_cp95"]
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    audit.check(len(rows) == 1 and set(rows[0]) == set(expected), "production_summary_csv_coverage", "all scalar fields")
    if len(rows) == 1:
        for key, wanted in expected.items():
            observed = rows[0].get(key)
            if isinstance(wanted, str):
                same = observed == wanted
            else:
                try:
                    same = math.isclose(float(observed), wanted, rel_tol=2e-10, abs_tol=2e-12)
                except (TypeError, ValueError):
                    same = False
            audit.check(same, "production_summary_csv_value", key)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--baseline-protocol", type=Path)
    parser.add_argument("--baseline-raw", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    original_protocol_path = args.baseline_protocol or root/"original"/"Protocol_JointRules.json"
    original_raw_path = args.baseline_raw or root/"original"/"raw"/"J1_n25600.csv.gz"
    protocol_path = root/"Protocol_J1_Replication.json"
    base = root/"results"/"main"
    raw_path = base/"raw"/"J1_n25600.csv.gz"
    detail_path = base/"details"/"J1_n25600.jsonl.gz"
    summary_path = base/"summary.json"
    csv_summary_path = base/"summary.csv"
    rule_count_path = base/"rule_counts.json"
    audit = Audit()
    for path in (Path(__file__).resolve(), original_protocol_path, original_raw_path, protocol_path, raw_path, detail_path, summary_path, csv_summary_path, rule_count_path):
        audit.hash(path, root)
    protocol = json.loads(protocol_path.read_text())
    original_protocol = json.loads(original_protocol_path.read_text())
    audit.check(protocol["fixed_design"] is True and protocol["sizes"] == [25600]
                and protocol["replicates_per_cell"] == protocol["total_draws"] == 20000
                and list(protocol["families"]) == ["J1"], "prespecified_scope", "one fixed cell")
    audit.check(protocol["master_seed"] == 2026101601 and protocol["families"]["J1"]["id"] == 3101,
                "prespecified_seed_namespace", "new master and family IDs")
    audit.check(protocol["master_seed"] != original_protocol["master_seed"], "nonoverlapping_seed_tuples", "distinct master")
    audit.check(protocol["component_caps"] == original_protocol["component_caps"], "unchanged_solver_caps", "same cap values")
    law = derive_law(protocol, original_protocol, audit)
    original_rows, replication_rows = read_rows(original_raw_path), read_rows(raw_path)
    check_rows(original_rows, original_protocol, audit, baseline=True)
    check_rows(replication_rows, protocol, audit)
    detail_accounting(detail_path, replication_rows, audit)
    original, replication = summarize(original_rows, "original joint-rule study"), summarize(replication_rows, "replication")
    audit.check(original["mean"]["integer_sums_lower"]["sum"] == 848
                and original["m2"]["integer_sums_lower"] == {"sum": 1162, "sum_of_squared_variables": 3742}
                and original["P0"]["definite_zero_count"] == 1293 and original["unresolved"] == 0,
                "fresh_original_baseline", "integer sums and zeros")
    production = json.loads(summary_path.read_text())
    audit.check(isinstance(production, list) and len(production) == 1, "production_cell_scope", "single cell summary")
    if isinstance(production, list) and production:
        compare_production(production[0], replication, law, audit)
    compare_summary_csv(csv_summary_path, replication, law, audit)
    diagnostic = rule_diagnostics(replication_rows)
    expected_rule_export = {"order": "corner_code0..7", "draws": 20000, "rows_per_degree": 256000000,
                            "counts": {str(item["degree"]): item["all_corner_counts"] for item in diagnostic["per_degree"]}}
    audit.check(json.loads(rule_count_path.read_text()) == expected_rule_export,
                "production_rule_counts_export", "both exact class totals and all eight categories")
    # Deliberately perturb only a diagnostic count, in memory. The real raw
    # records are never changed, and this record is not another sampled draw.
    perturbed = [dict(row) for row in replication_rows]
    changed_counts = json.loads(perturbed[0]["rule_counts_d1"])
    changed_counts[0] -= 1
    changed_counts[5] += 1
    perturbed[0]["rule_counts_d1"] = json.dumps(changed_counts)
    invariant_cost = summarize(perturbed, "replication") == replication
    invariant_reference = derive_law(protocol, original_protocol, audit) == law
    diagnostic_changed = rule_diagnostics(perturbed) != diagnostic
    audit.check(invariant_cost and invariant_reference and diagnostic_changed,
                "no_refit_mutation_check", "one ZERO count moved to COPY1 in an in-memory diagnostic-only record")
    report = {"status": "FAIL" if audit.failures else "PASS", "schema_version": 1,
              "scope": "Independent raw/detail accounting, exact-integer aggregation, fixed-law derivation and prespecified sampling diagnostics. Network regeneration and witness dynamics are verified separately.",
              "method": {"normal_quantile": Z95, "moment_intervals": "Pointwise approximate normal envelopes for finite-n endpoint moments; exact integer sums precede division. Computational bounds and Monte Carlo uncertainty are distinct.",
                         "proportion_intervals": "Pointwise equal-tailed 95% Clopper-Pearson intervals; outer envelopes if exact cost is unresolved.",
                         "reference": "Fixed limiting Poisson reference, never fitted. Coverage is not a finite-size approximation-error bound."},
              "fixed_reference": law, "original": original, "replication": replication,
              "no_refit_mutation_check": {"cost_summary_unchanged": invariant_cost,
                                           "fixed_reference_unchanged": invariant_reference,
                                           "diagnostics_changed": diagnostic_changed,
                                           "input_files_modified": False},
              "batch_difference": differences(original, replication), "rule_diagnostics": diagnostic,
              "checks": dict(sorted(audit.checks.items())), "failure_count": len(audit.failures),
              "failures": audit.failures, "source_sha256": dict(sorted(audit.hashes.items()))}
    output = args.output or root/"qa"/"statistical_result.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+"\n")
    flat_rows = []
    for summary in (original, replication):
        for variable, reference in (("mean", law["mu"]), ("m2", law["reference_m2"]), ("P0", law["reference_p0"])):
            item = summary[variable]
            interval = item["cp95_envelope"] if variable == "P0" else item["normal95_envelope"]
            flat_rows.append({"sample": summary["label"], "draws": summary["draws"], "n": summary["n"],
                              "quantity": variable, "estimate_lower": item["lower"], "estimate_upper": item["upper"],
                              "se_lower": item.get("se_lower"), "se_upper": item.get("se_upper"),
                              "interval95_lower": interval[0], "interval95_upper": interval[1],
                              "interval_type": "Clopper-Pearson" if variable == "P0" else "approximate normal",
                              "limiting_reference": reference, "unresolved": summary["unresolved"]})
    write_csv(base/"statistical_summary.csv", flat_rows)
    write_csv(base/"rule_frequency_summary.csv", diagnostic["rows"])
    print(json.dumps({"status": report["status"], "draws": replication["draws"], "failures": audit.failures[:10],
                      "mean": replication["mean"], "m2": replication["m2"], "P0": replication["P0"],
                      "joint_rule_diagnostic_p": diagnostic["joint_diagnostic_p_value"]}, indent=2))
    return bool(audit.failures)


if __name__ == "__main__":
    raise SystemExit(main())
