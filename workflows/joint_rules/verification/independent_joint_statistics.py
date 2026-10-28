#!/usr/bin/env python3
"""Independent complete-sample statistical audit for the focused JointRules study.

Run from any directory: python verification/independent_joint_statistics.py
No production modules are imported. Law parameters are derived with exact
fractions from the declared corner-code weights; all empirical moments start
with integer sums. SciPy is used only as the beta-quantile numerical oracle.
This verifies aggregation, accounting and provenance, not network regeneration
or the validity of a dynamic control witness (checked by separate validators).
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist
import sys

from scipy.special import betaincinv


CP_CACHE = {}


class Audit:
    def __init__(self):
        self.checks = Counter()
        self.failures = []
        self.sources = {}
        self.largest_float_error = 0.0

    def check(self, condition, category, label):
        self.checks[category] += 1
        if not condition:
            self.failures.append(label)

    def compare(self, expected, actual, label):
        if isinstance(expected, dict):
            if isinstance(actual, str):
                try:
                    actual = json.loads(actual)
                except ValueError:
                    actual = None
            self.check(isinstance(actual, dict), "export_field_types", label)
            if not isinstance(actual, dict):
                return
            self.check(set(expected) == set(actual), "export_field_coverage", label)
            for key, value in expected.items():
                self.compare(value, actual.get(key), label+"/"+key)
        elif isinstance(expected, list):
            self.check(isinstance(actual, list) and len(actual) == len(expected), "export_array_coverage", label)
            if isinstance(actual, list):
                for i, (wanted, observed) in enumerate(zip(expected, actual)):
                    self.compare(wanted, observed, f"{label}/{i}")
        elif isinstance(expected, int) and not isinstance(expected, bool):
            try:
                same = Decimal(str(actual)) == Decimal(expected)
            except (InvalidOperation, TypeError):
                same = False
            self.check(same, "export_exact_integer_comparisons", label)
        elif isinstance(expected, float):
            try:
                value = float(actual)
                difference = abs(expected-value)
                self.largest_float_error = max(self.largest_float_error, difference)
                same = math.isfinite(value) and math.isclose(expected, value, rel_tol=2e-10, abs_tol=2e-12)
            except (ValueError, TypeError):
                same = False
            self.check(same, "export_float_comparisons", label)
        elif expected is None:
            self.check(actual in (None, ""), "export_other_comparisons", label)
        else:
            self.check(actual == expected, "export_other_comparisons", label)

    def hash(self, path, root):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        self.sources[str(path.relative_to(root))] = digest
        return digest


def fraction_pair(value):
    return [value.numerator, value.denominator]


def read_csv(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream))


def cp95(k, n):
    if (k, n) not in CP_CACHE:
        CP_CACHE[k, n] = [0.0 if k == 0 else float(betaincinv(k, n-k+1, .025)),
                          1.0 if k == n else float(betaincinv(k+1, n-k, .975))]
    return CP_CACHE[k, n]


def binomial_probability_interval(n, p, first, last):
    logp, logq = math.log(p), math.log1p(-p)
    logterm = (math.lgamma(n+1)-math.lgamma(first+1)-math.lgamma(n-first+1)
               + first*logp+(n-first)*logq)
    terms = [math.exp(logterm)]
    for k in range(first, last):
        logterm += math.log(n-k)-math.log(k+1)+logp-logq
        terms.append(math.exp(logterm))
    return math.fsum(terms)


def check_cp_quantiles(audit):
    worst = 0.0
    for (k, n), (low, high) in sorted(CP_CACHE.items()):
        if k:
            error = abs(binomial_probability_interval(n, low, k, n)-.025)
            worst = max(worst, error)
            audit.check(error <= 2e-10, "CP_binomial_tail_identity", f"CP lower/{k}/{n}")
        if k < n:
            error = abs(binomial_probability_interval(n, high, 0, k)-.025)
            worst = max(worst, error)
            audit.check(error <= 2e-10, "CP_binomial_tail_identity", f"CP upper/{k}/{n}")
    return {"distinct_count_intervals": len(CP_CACHE), "maximum_absolute_tail_error": worst,
            "absolute_tail_tolerance": 2e-10}


def parameters(protocol, audit):
    results = {}
    for name, family in protocol["families"].items():
        lam = rho = mass = mean_degree = Fraction(0)
        type_records = []
        for index, type_info in enumerate(family["types"]):
            tag = f"law/{name}/type/{index}"
            fraction = Fraction(*type_info["fraction"])
            weights = type_info["weights"]
            denominator = type_info["probability_denominator"]
            degree = type_info["degree"]
            audit.check(len(weights) == 8 and all(isinstance(w, int) and w >= 0 for w in weights)
                        and sum(weights) == denominator and denominator > 0, "exact_corner_probability_law", tag)
            probability = [Fraction(w, denominator) for w in weights]
            q1 = sum((p for code, p in enumerate(probability) if code & 1), Fraction(0))
            q2 = sum((p for code, p in enumerate(probability) if code & 2), Fraction(0))
            q12 = sum((p for code, p in enumerate(probability) if code & 3 == 3), Fraction(0))
            nonzero = 1-probability[0]
            derived = {"q1": fraction_pair(q1), "q2": fraction_pair(q2),
                       "q12": fraction_pair(q12), "nonzero_row_probability": fraction_pair(nonzero)}
            audit.compare(derived, type_info["derived_rational"], tag+"/derived_rational")
            audit.check(fraction > 0 and isinstance(degree, int) and degree >= 0, "valid_type_fraction_degree", tag)
            for size in protocol["sizes"]:
                audit.check((fraction*size).denominator == 1, "integer_type_block_sizes", tag+f"/n{size}")
            mass += fraction
            mean_degree += fraction*degree
            lam += fraction*degree*(q1+q2)/2
            rho += fraction*degree*nonzero/2
            # Directly average support over four independent free-slot patterns.
            support = [Fraction(0), q1, q2, nonzero]
            type_records.append({"degree": degree, "fraction": fraction_pair(fraction),
                                 "corner_probabilities": [fraction_pair(p) for p in probability],
                                 "derived_rational": derived,
                                 "support_probabilities_masks_0_1_2_3": [fraction_pair(p) for p in support]})
        audit.check(mass == 1 and mean_degree == 2, "law_normalization_and_stub_balance", name)
        contraction = max(lam, rho)
        audit.check(0 <= lam < 1 and 0 <= rho < 1, "strict_contraction", name)
        expected = {"lambda": fraction_pair(lam), "rho": fraction_pair(rho),
                    "contraction_coefficient": fraction_pair(contraction)}
        audit.compare(expected, family["parameters_rational"], name+"/parameters_rational")
        mu = -math.log1p(-float(lam))
        audit.compare(mu, family["reference_mu"], name+"/reference_mu")
        psi_grid = []
        for numerator in range(17):
            u = Fraction(numerator, 16)
            direct = Fraction(0)
            for info in type_records:
                q1 = Fraction(*info["derived_rational"]["q1"])
                q2 = Fraction(*info["derived_rational"]["q2"])
                r = Fraction(*info["derived_rational"]["nonzero_row_probability"])
                direct += Fraction(*info["fraction"])*info["degree"]*((q1+q2)*u*(1-u)+r*u*u)/2
            polynomial = lam*u+(rho-lam)*u*u
            audit.check(direct == polynomial and polynomial <= contraction*u, "exact_support_map_identity", f"{name}/u={u}")
            psi_grid.append({"u": fraction_pair(u), "psi": fraction_pair(polynomial)})
        results[name] = {"parameters_rational": expected, "types": type_records,
                         "psi_coefficients_u_u2": [fraction_pair(lam), fraction_pair(rho-lam)],
                         "psi_rational_grid": psi_grid, "mu": mu, "reference_m2": mu+mu*mu,
                         "reference_p0": float(1-lam)}
    return results


def moment(values, power):
    count = len(values)
    total = sum(x**power for x in values)
    squares = sum(x**(2*power) for x in values)
    mean = total/count
    se = math.sqrt((count*squares-total*total)/(count*count*(count-1)))
    return mean, se, total, squares


def summarize(rows, law, z):
    count = len(rows)
    lower = [int(row["cost_lower"]) for row in rows]
    upper = [int(row["cost_upper"]) for row in rows]
    gfail = sum(int(row["G"]) == 0 for row in rows)
    gpfail = sum(int(row["G_prime"]) == 0 for row in rows)
    definite, possible = upper.count(0), lower.count(0)
    output = {"family": rows[0]["family"], "n": int(rows[0]["n"]), "draws": count,
              "lambda": float(Fraction(*law["parameters_rational"]["lambda"])),
              "rho": float(Fraction(*law["parameters_rational"]["rho"])),
              "contraction_coefficient": float(Fraction(*law["parameters_rational"]["contraction_coefficient"])),
              "mu": law["mu"], "reference_m2": law["reference_m2"], "reference_p0": law["reference_p0"],
              "unresolved": sum(lo != hi for lo, hi in zip(lower, upper)),
              "G_failures": gfail, "G_prime_failures": gpfail,
              "G_failures_cp95": cp95(gfail, count), "G_prime_failures_cp95": cp95(gpfail, count),
              "P0": {"definite_zero_count": definite, "possible_zero_count": possible,
                     "lower": definite/count, "upper": possible/count,
                     "cp95_envelope": [cp95(definite, count)[0], cp95(possible, count)[1]]}}
    integer_sums = {}
    for label, power in (("mean", 1), ("m2", 2)):
        ml, sl, suml, squaresl = moment(lower, power)
        mu, su, sumu, squaresu = moment(upper, power)
        output[label] = {"lower": ml, "upper": mu, "se_lower": sl, "se_upper": su,
                         "normal95_envelope": [max(0.0, ml-z*sl), mu+z*su]}
        integer_sums[label] = {"lower_sum": suml, "lower_sum_squared_variables": squaresl,
                               "upper_sum": sumu, "upper_sum_squared_variables": squaresu}
    stops = Counter()
    for row in rows:
        stops.update(json.loads(row["stop_summary"]))
    reasons = dict(sorted(Counter(row["exact_reason"] for row in rows).items()))
    output.update({"max_residual": max(int(row["residual_size"]) for row in rows),
                   "max_complex": max(int(row["largest_complex"]) for row in rows),
                   "total_solver_work": sum(int(row["solver_work"]) for row in rows),
                   "stop_reason_counts": dict(sorted(stops.items())), "exact_reason_counts": reasons})
    extras = {"integer_moment_sums": integer_sums,
              "G_prime_rescued_count": sum(int(row["G_prime"]) == 1 and int(row["G"]) == 0 for row in rows),
              "exact_reason_counts": reasons,
              "stop_reason_component_counts": dict(sorted(stops.items())),
              "residual_size_histogram": dict(sorted(Counter(int(row["residual_size"]) for row in rows).items())),
              "largest_complex_histogram": dict(sorted(Counter(int(row["largest_complex"]) for row in rows).items())),
              "cost_interval_width_histogram": dict(sorted(Counter(hi-lo for lo, hi in zip(lower, upper)).items())),
              "solver_work_total": sum(int(row["solver_work"]) for row in rows)}
    return output, extras


def check_rows(rows, family, size, replicates, audit):
    tag = f"{family}_n{size}"
    audit.check(len(rows) == replicates, "complete_cell_row_count", tag)
    audit.check(sorted(int(row["replicate"]) for row in rows) == list(range(replicates)),
                "complete_unique_replicate_ids", tag)
    for row in rows:
        label = tag+"/"+row["replicate"]
        audit.check(row["family"] == family and int(row["n"]) == size, "row_cell_identity", label)
        lo, hi, residual = (int(row[key]) for key in ("cost_lower", "cost_upper", "residual_size"))
        audit.check(0 <= lo <= hi <= residual <= size, "valid_cost_interval", label)
        g, gp = int(row["G"]), int(row["G_prime"])
        audit.check(g in (0, 1) and gp in (0, 1) and (not g or gp), "nested_certificates", label)
        audit.check(g == (int(row["double_rows"]) == 0), "G_double_row_equivalence", label)
        audit.check(not gp or lo == hi == int(row["certificate_cycle_count"]), "certificate_cost_identity", label)
        audit.check(gp or row["certificate_cycle_count"] == "", "uncertified_cycle_count_representation", label)
        pins = [int(pin) for pin in row["pins"].split(";") if pin]
        audit.check(len(pins) == len(set(pins)) == hi and all(0 <= pin < size for pin in pins),
                    "witness_cardinality_and_labels", label)
        digest = row["network_sha256"]
        audit.check(len(digest) == 64 and all(c in "0123456789abcdef" for c in digest), "network_hash_format", label)
        stops = json.loads(row["stop_summary"])
        audit.check(isinstance(stops, dict) and all(isinstance(value, int) and value >= 0 for value in stops.values()),
                    "stop_count_format", label)
        expected_reason = ("original_G" if g else "weaker_G_prime" if gp else
                           "SCC_exact" if lo == hi else "bounded_unresolved")
        audit.compare(expected_reason, row["exact_reason"], "row_exact_reason/"+label)
        audit.check(0 <= int(row["cyclic_components"]) <= residual and 0 <= int(row["largest_complex"]) <= residual
                    and int(row["solver_work"]) >= 0, "SCC_and_work_ranges", label)


def compare_csv(path, expected, audit, root):
    audit.hash(path, root)
    observed = read_csv(path)
    lookup = {(str(row["family"]), int(row["n"])): row for row in expected}
    observed_keys = [(row["family"], int(row["n"])) for row in observed]
    audit.check(set(observed_keys) == set(lookup) and len(observed_keys) == len(lookup), "summary_csv_cell_coverage", path.name)
    for row in observed:
        key = row["family"], int(row["n"])
        if key in lookup:
            audit.compare(lookup[key], row, f"{path.name}/{key}")


def flattened(summary):
    result = {key: value for key, value in summary.items() if not isinstance(value, (dict, list))}
    for variable in ("mean", "m2"):
        for field in ("lower", "upper", "se_lower", "se_upper"):
            result[f"{variable}_{field}"] = summary[variable][field]
        result[f"{variable}_normal95_lo"], result[f"{variable}_normal95_hi"] = summary[variable]["normal95_envelope"]
    for field in ("definite_zero_count", "possible_zero_count", "lower", "upper"):
        result["P0_"+field] = summary["P0"][field]
    result["P0_cp95_lo"], result["P0_cp95_hi"] = summary["P0"]["cp95_envelope"]
    for field in ("G_failures", "G_prime_failures"):
        result[field+"_cp95_lo"], result[field+"_cp95_hi"] = summary[field+"_cp95"]
    return result


def check_detail_stream(path, raw_rows, audit):
    """Transport/accounting check only; no dynamic claim follows from this."""
    count = 0
    with gzip.open(path, "rt") as stream:
        for count, line in enumerate(stream, 1):
            detail = json.loads(line)
            if count > len(raw_rows):
                audit.check(False, "detail_row_identity", f"{path.name}/extra/{count}")
                continue
            raw = raw_rows[count-1]
            identity = f"{raw['family']}/{raw['n']}/{raw['replicate']}"
            same = (detail["family"] == raw["family"] and detail["n"] == int(raw["n"])
                    and detail["replicate"] == int(raw["replicate"])
                    and detail["network_sha256"] == raw["network_sha256"])
            audit.check(same, "detail_row_identity", identity)
            pins = [int(v) for v in raw["pins"].split(";") if v]
            same = (detail["lower"] == int(raw["cost_lower"]) and detail["upper"] == int(raw["cost_upper"])
                    and detail["G"] == bool(int(raw["G"])) and detail["G_prime"] == bool(int(raw["G_prime"]))
                    and detail["exact_reason"] == raw["exact_reason"] and detail["pins"] == pins
                    and len(detail["residual_rows"]) == int(raw["residual_size"]))
            audit.check(same, "detail_raw_summary_consistency", identity)
            expected_cycle_count = int(raw["certificate_cycle_count"]) if int(raw["G_prime"]) else None
            audit.check(detail["certificate_cycle_count"] == expected_cycle_count, "detail_raw_cycle_count_consistency", identity)
    audit.check(count == len(raw_rows), "detail_stream_record_count", path.name)
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output or root/"qa"/"independent_joint_statistics.json"
    validation = root/"validation"
    base = validation/"results"/"main"
    protocol_path = validation/"Protocol_JointRules.json"
    audit = Audit()
    protocol = json.loads(protocol_path.read_text())
    protocol_hash = audit.hash(protocol_path, root)
    audit.hash(Path(__file__).resolve(), root)
    audit.check(protocol["fixed_design"] is True, "fixed_design", "fixed_design")
    audit.check(len(protocol["families"]) == 3 and len(protocol["sizes"]) == 3
                and protocol["replicates_per_cell"] == 2000 and protocol["total_draws"] == 18000,
                "focused_design_scope", "3 laws x 3 sizes x 2000 draws")
    family_ids = [info["id"] for info in protocol["families"].values()]
    audit.check(len(family_ids) == len(set(family_ids)), "unique_seed_family_ids", "family IDs")
    law_parameters = parameters(protocol, audit)
    vendor_sources = {"vendor/engine.py", "vendor/engine_core.cpp", "vendor/scc_control.py"}
    audit.check(set(protocol["source_sha256"]) == vendor_sources, "source_hash_key_coverage", "protocol source keys")
    for relative, recorded_hash in protocol["source_sha256"].items():
        audit.compare(audit.hash(validation/relative, root), recorded_hash, "protocol/source/"+relative)
    manifest_path = base/"manifest.json"
    if not manifest_path.is_file():
        parser.error("The completed 18,000-draw focused main replay is required.")
    audit.hash(manifest_path, root)
    manifest = json.loads(manifest_path.read_text())
    audit.compare("completed", manifest["status"], "manifest/status")
    audit.compare(protocol["total_draws"], manifest["total_draws"], "manifest/total_draws")
    audit.compare(protocol_hash, manifest["protocol_sha256"], "manifest/protocol_sha256")
    gate_path = root/"qa"/"independent_exact_checks.json"
    audit.compare(audit.hash(gate_path, root), manifest["independent_gate_sha256"], "manifest/independent_gate_sha256")
    gate = json.loads(gate_path.read_text())
    audit.compare("PASS", gate["status"], "independent_gate/status")
    audit.compare(protocol_hash, gate["protocol_sha256"], "independent_gate/protocol_sha256")
    audit.check(set(gate["source_hashes"]) == {"independent_exact_checks.py", "independent_network_oracle.py"},
                "source_hash_key_coverage", "independent gate source keys")
    for filename, recorded_hash in gate["source_hashes"].items():
        audit.compare(audit.hash(root/"verification"/filename, root), recorded_hash, "independent_gate/source/"+filename)
    audit.check(set(manifest["source_sha256"]) == vendor_sources | {"src/joint_validation.py"},
                "source_hash_key_coverage", "manifest source keys")
    for relative, recorded_hash in manifest["source_sha256"].items():
        audit.compare(audit.hash(validation/relative, root), recorded_hash, "manifest/source/"+relative)
    cells = [(name, size) for name in sorted(protocol["families"]) for size in protocol["sizes"]]
    replicates = protocol["replicates_per_cell"]
    raw_names = {f"{name}_n{size}.csv.gz" for name, size in cells}
    detail_names = {f"{name}_n{size}.jsonl.gz" for name, size in cells}
    audit.check({path.name for path in (base/"raw").glob("*.csv.gz")} == raw_names, "complete_raw_cell_files", "raw cells")
    audit.check({path.name for path in (base/"details").glob("*.jsonl.gz")} == detail_names, "complete_detail_cell_files", "detail cells")
    z = NormalDist().inv_cdf(.975)
    expected_summaries = []
    extra_statistics = []
    expected_integrity = []
    seed_identity = hashlib.sha256()
    total_raw_rows = total_detail_rows = 0
    for name, size in cells:
        raw_path = base/"raw"/f"{name}_n{size}.csv.gz"
        detail_path = base/"details"/f"{name}_n{size}.jsonl.gz"
        raw_hash = audit.hash(raw_path, root)
        rows = read_csv(raw_path)
        check_rows(rows, name, size, replicates, audit)
        total_raw_rows += len(rows)
        detail_hash = audit.hash(detail_path, root)
        total_detail_rows += check_detail_stream(detail_path, rows, audit)
        # Check hashes again after complete decompression to detect changing inputs.
        audit.compare(raw_hash, hashlib.sha256(raw_path.read_bytes()).hexdigest(), "stable_raw_stream/"+raw_path.name)
        audit.compare(detail_hash, hashlib.sha256(detail_path.read_bytes()).hexdigest(), "stable_detail_stream/"+detail_path.name)
        for row in sorted(rows, key=lambda r: int(r["replicate"])):
            seed_identity.update(f"{protocol['master_seed']}:{protocol['families'][name]['id']}:{size}:{row['replicate']}\n".encode())
        summary, extras = summarize(rows, law_parameters[name], z)
        summary["stream_sha256"] = {"raw": raw_hash, "details": detail_hash}
        expected_summaries.append(summary)
        extra_statistics.append({"family": name, "n": size, **extras})
        for kind, path, digest in (("raw", raw_path, raw_hash), ("details", detail_path, detail_hash)):
            expected_integrity.append({"path": str(path.relative_to(base)), "records": len(rows), "sha256": digest})
    summary_path = base/"summary.json"
    audit.hash(summary_path, root)
    audit.compare(expected_summaries, json.loads(summary_path.read_text()), "summary.json")
    csv_rows = [flattened(summary) for summary in expected_summaries]
    compare_csv(base/"summary.csv", csv_rows, audit, root)
    integrity_path = base/"io_integrity.json"
    audit.hash(integrity_path, root)
    audit.compare({"status": "PASS", "streams": expected_integrity}, json.loads(integrity_path.read_text()), "io_integrity.json")
    total_unresolved = sum(summary["unresolved"] for summary in expected_summaries)
    audit.compare(total_unresolved, manifest["unresolved"], "manifest/unresolved")
    audit.check(total_raw_rows == total_detail_rows == protocol["total_draws"], "complete_sample_accounting", "all raw/detail rows")
    # Include both exact binomial endpoint cases even if a particular dataset has none.
    cp95(0, replicates)
    cp95(replicates, replicates)
    quantile_validation = check_cp_quantiles(audit)
    report = {
        "schema_version": 1, "status": "FAIL" if audit.failures else "PASS",
        "scope": "Independent rational-law derivation, complete-sample arithmetic, raw/detail accounting and provenance; no production imports, network generation, state-graph solver, or dynamic witness validation.",
        "coverage": {"families": len(protocol["families"]), "cells": len(cells), "unique_draws": total_raw_rows,
                     "raw_rows_read": total_raw_rows, "detail_rows_read": total_detail_rows,
                     "summary_JSON_cells": len(expected_summaries), "summary_CSV_cells": len(csv_rows),
                     "summary_CSV_fields_per_cell": len(csv_rows[0])},
        "method": {"moments": "Exact integer first/second/fourth sums before floating division; unbiased sample variance divided by sample count for each endpoint variable.",
                   "normal_quantile": z,
                   "sampling_envelopes": "Pointwise approximate normal outer envelopes for raw moments; pointwise Clopper-Pearson outer envelope for zero-cost probability. Computational intervals are separate from sampling uncertainty; no imputation or exclusions.",
                   "CP_numerics": "scipy.special.betaincinv, independently checked by direct binomial log-pmf summation at both tails; boundary counts 0 and N included.",
                   "numeric_comparison_tolerance": {"absolute": 2e-12, "relative": 2e-10},
                   "seed_identity": "Complete protocol master_seed/family_id/n/replicate tuples; this is record identity, not RNG replay.",
                   "certificate_semantics": "Nominal residual slots, retaining duplicate parents as separate certificate slots."},
        "law_parameters": law_parameters,
        "totals": {"draws": total_raw_rows, "exact": total_raw_rows-total_unresolved, "unresolved": total_unresolved,
                   "G_failures": sum(summary["G_failures"] for summary in expected_summaries),
                   "G_prime_failures": sum(summary["G_prime_failures"] for summary in expected_summaries),
                   "G_prime_rescued_count": sum(item["G_prime_rescued_count"] for item in extra_statistics)},
        "expected_summary": expected_summaries, "additional_statistics": extra_statistics,
        "seed_identity_sha256": seed_identity.hexdigest(), "source_sha256": dict(sorted(audit.sources.items())),
        "checks": dict(sorted(audit.checks.items())), "CP_quantile_validation": quantile_validation,
        "maximum_absolute_export_float_error": audit.largest_float_error,
        "failure_count": len(audit.failures), "failures": audit.failures,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, sort_keys=True, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"status": report["status"], "coverage": report["coverage"], "totals": report["totals"],
                      "failure_count": len(audit.failures), "first_failures": audit.failures[:12]}, indent=2))
    return 1 if audit.failures else 0


if __name__ == "__main__":
    sys.exit(main())
