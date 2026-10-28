#!/usr/bin/env python3
"""Independently aggregate all Baseline/Bernoulli raw records without importing study code.

Run from any directory:
  python verification/independent_statistics.py
SciPy supplies beta quantiles and a separate external oracle for exported
confidence limits and chi-square probabilities. Integer histograms and sums
determine every empirical statistic. This checks record
identity, arithmetic, interval nesting, and exported baseline aggregates; it does
not regenerate a network or independently certify a control witness.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys

from scipy.special import betaincinv
from scipy.stats import beta as scipy_beta, chi2 as scipy_chi2


CP_ENDPOINTS = {}


class Audit:
    def __init__(self):
        self.checks = Counter()
        self.failures = []
        self.source_hashes = {}
        self.scipy_errors = Counter()

    def check(self, condition, category, label):
        self.checks[category] += 1
        if not condition:
            self.failures.append(label)

    def compare(self, expected, actual, label):
        if isinstance(expected, int) and not isinstance(expected, bool):
            try:
                good = Decimal(str(actual)) == Decimal(expected)
            except (InvalidOperation, TypeError):
                good = False
        elif isinstance(expected, float):
            try:
                value = float(actual)
                good = math.isfinite(value) and math.isclose(
                    expected, value, rel_tol=2e-10, abs_tol=2e-12)
            except (ValueError, TypeError):
                good = False
        elif expected is None:
            good = actual in (None, "")
        elif isinstance(expected, dict) and isinstance(actual, str):
            try:
                good = expected == json.loads(actual)
            except ValueError:
                good = False
        else:
            good = expected == actual
        self.check(good, "exported_statistic_comparisons", label)

    def hash(self, path, root):
        self.source_hashes[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()


def compare_scipy_value(expected, actual, category, label, audit):
    observed = float(actual)
    error = abs(expected-observed)
    audit.scipy_errors[category] = max(audit.scipy_errors[category], error)
    audit.check(math.isfinite(observed) and math.isclose(expected, observed, rel_tol=2e-10, abs_tol=2e-12),
                category, label)


def compare_scipy_cp(k, n, endpoint, side, label, audit):
    if side == "lower":
        reference = 0.0 if k == 0 else float(scipy_beta.ppf(0.025, k, n-k+1))
    else:
        reference = 1.0 if k == n else float(scipy_beta.ppf(0.975, k+1, n-k))
    compare_scipy_value(reference, endpoint, "SciPy_beta_ppf_export_endpoint", label, audit)


def compare_scipy_summary(summary, label, audit):
    n = int(summary["M"])
    for event in ("g", "g_prime"):
        k = int(summary[event+"_failures"])
        compare_scipy_cp(k, n, summary[event+"_cp_lo"], "lower", label+"/"+event+"/lower", audit)
        compare_scipy_cp(k, n, summary[event+"_cp_hi"], "upper", label+"/"+event+"/upper", audit)
    compare_scipy_cp(int(summary["p0_definite_count"]), n, summary["p0_cp_lo"], "lower", label+"/p0/lower", audit)
    compare_scipy_cp(int(summary["p0_possible_count"]), n, summary["p0_cp_hi"], "upper", label+"/p0/upper", audit)
    for side, stat_side in (("lower", "upper"), ("upper", "lower")):
        reference = float(scipy_chi2.sf(float(summary["chi_square_statistic_"+stat_side]), int(summary["chi_square_df"])))
        compare_scipy_value(reference, summary["chi_square_p_"+side], "SciPy_chi2_sf_export_bound", label+"/chi_square/"+side, audit)


def read_csv(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream))


def cp95(k, n):
    if (k, n) not in CP_ENDPOINTS:
        CP_ENDPOINTS[k, n] = [0.0 if k == 0 else float(betaincinv(k, n-k+1, 0.025)),
                              1.0 if k == n else float(betaincinv(k+1, n-k, 0.975))]
    return CP_ENDPOINTS[k, n]


def binomial_interval_probability(n, p, first, last):
    """Direct log-pmf summation, independent of incomplete beta functions."""
    logp, logq = math.log(p), math.log1p(-p)
    term = (math.lgamma(n+1)-math.lgamma(first+1)-math.lgamma(n-first+1)
            + first*logp + (n-first)*logq)
    values = [math.exp(term)]
    for j in range(first, last):
        term += math.log(n-j)-math.log(j+1)+logp-logq
        values.append(math.exp(term))
    return math.fsum(values)


def verify_cp_identities(audit):
    maximum_error = 0.0
    for (k, n), (lo, hi) in sorted(CP_ENDPOINTS.items()):
        if k:
            error = abs(binomial_interval_probability(n, lo, k, n)-0.025)
            maximum_error = max(maximum_error, error)
            audit.check(error <= 2e-10, "CP_binomial_tail_identity", f"lower/{k}/{n}/{error}")
        if k < n:
            error = abs(binomial_interval_probability(n, hi, 0, k)-0.025)
            maximum_error = max(maximum_error, error)
            audit.check(error <= 2e-10, "CP_binomial_tail_identity", f"upper/{k}/{n}/{error}")
    return {"unique_intervals": len(CP_ENDPOINTS), "maximum_absolute_tail_probability_error": maximum_error,
            "absolute_tolerance": 2e-10}


def chi_square_sf(statistic, df):
    """Gamma recurrence for integer/half-integer df/2, using math.erfc only."""
    if statistic <= 0:
        return 1.0
    x = statistic/2
    if df % 2:
        value, shape = math.erfc(math.sqrt(x)), 0.5
    else:
        value, shape = math.exp(-x), 1.0
    while shape < df/2:
        value += math.exp(shape*math.log(x)-x-math.lgamma(shape+1))
        shape += 1
    return min(1.0, max(0.0, value))


def sample_mean_se(histogram, power=1):
    """Exact integer sums before the final division and square root."""
    n = sum(histogram.values())
    first = sum(count * value**power for value, count in histogram.items())
    second = sum(count * value**(2*power) for value, count in histogram.items())
    variance_of_mean = (n*second - first*first) / (n*n*(n-1)) if n > 1 else 0.0
    return first/n, math.sqrt(variance_of_mean), first, second


def poisson_sequence(mean, maximum):
    probabilities = [math.exp(-mean)]
    for k in range(1, maximum+1):
        probabilities.append(probabilities[-1]*mean/k)
    cdf = [min(1.0, math.fsum(probabilities[:k+1])) for k in range(maximum+1)]
    return probabilities, cdf


def aggregate(rows, mixture, size, lam, stage):
    n = len(rows)
    lower = Counter(int(r["cost_lower"]) for r in rows)
    upper = Counter(int(r["cost_upper"]) for r in rows)
    residual = Counter(int(r["residual_size"]) for r in rows)
    exact = Counter(int(r["cost_lower"]) for r in rows if r["cost_lower"] == r["cost_upper"])
    widths = Counter(int(r["cost_upper"])-int(r["cost_lower"]) for r in rows)
    maximum = max(upper)
    mean = -math.log1p(-lam)
    poisson_pmf, poisson_cdf = poisson_sequence(mean, max(maximum, 11))
    distribution = []
    cumulative_l = cumulative_u = 0
    discrepancy_l = discrepancy_u = 0.0
    for k in range(maximum+1):
        cumulative_l += upper[k]  # definitely at or below k
        cumulative_u += lower[k]  # possibly at or below k
        possible = sum(count for value, count in lower.items() if value <= k) - sum(
            count for value, count in upper.items() if value < k)
        fl, fu = cumulative_l/n, cumulative_u/n
        target = poisson_cdf[k]
        discrepancy_l = max(discrepancy_l, fl-target, target-fu)
        discrepancy_u = max(discrepancy_u, abs(fl-target), abs(fu-target))
        distribution.append({"k": k, "pmf_lower": exact[k]/n,
                             "pmf_upper": possible/n,
                             "cdf_lower": fl, "cdf_upper": fu,
                             "poisson_pmf": poisson_pmf[k], "poisson_cdf": target})
    bins = []
    for k in range(12):
        if k < 11:
            definite = exact[k]
            possible = sum(count for value, count in lower.items() if value <= k) - sum(
                count for value, count in upper.items() if value < k)
            target = poisson_pmf[k]
        else:
            definite = sum(count for value, count in lower.items() if value >= 11)
            possible = sum(count for value, count in upper.items() if value >= 11)
            target = max(0.0, 1-poisson_cdf[10])
        bins.append({"bin": k, "pmf_lower": definite/n, "pmf_upper": possible/n,
                     "poisson_pmf": target})
    survival = []
    remaining = n
    for threshold in range(max(residual)+1):
        survival.append({"threshold": threshold, "residual_survival_ge": remaining/n})
        remaining -= residual[threshold]
    residual_mean, residual_se, _, _ = sample_mean_se(residual)
    certified = sum(int(r["G"] if stage == "bernoulli" else r["certified"]) for r in rows)
    failures = n-certified
    cp_l, cp_u = cp95(failures, n)
    work = [int(r["solver_work"]) for r in rows]
    resolved = sum(exact.values())
    summary = {
        "mixture": mixture, "n": size, "draws": n, "lambda_parameter": lam,
        "poisson_mean": mean, "poisson_raw_second_moment": mean+mean*mean,
        "certified_count": certified, "certificate_failure_count": failures,
        "certificate_failure_fraction": failures/n,
        "certificate_failure_cp95_lower": cp_l, "certificate_failure_cp95_upper": cp_u,
        "resolved_count": resolved, "unresolved_count": n-resolved,
        "resolved_fraction": resolved/n,
        "residual_mean": residual_mean, "residual_mean_mc_se": residual_se,
        "residual_max": max(residual), "double_rows_total": sum(int(r["double_rows"]) for r in rows),
        "solver_work_total": sum(work), "solver_work_max": max(work),
        "empirical_cdf_discrepancy_lower": discrepancy_l,
        "empirical_cdf_discrepancy_upper": discrepancy_u,
    }
    integer_sums = {}
    for power in (1, 2):
        for side, histogram in (("lower", lower), ("upper", upper)):
            value, se, sum1, sum2 = sample_mean_se(histogram, power)
            summary[f"raw_moment_{power}_{side}"] = value
            summary[f"endpoint_raw_moment_{power}_{side}_mc_se"] = se
            integer_sums[f"sum_{side}_power_{power}"] = sum1
            integer_sums[f"sum_{side}_power_{2*power}"] = sum2
        summary[f"resolved_raw_moment_{power}_mc_se"] = (
            summary[f"endpoint_raw_moment_{power}_lower_mc_se"] if resolved == n else None)
    if stage == "bernoulli":
        gp = sum(int(r["G_prime"]) for r in rows)
        gp_l, gp_u = cp95(n-gp, n)
        reason = Counter(r["exact_reason"] for r in rows)
        stops = Counter()
        for r in rows:
            stops.update(json.loads(r["stop_summary"]))
        summary.update({
            "G_count": certified, "G_failure_count": failures,
            "G_prime_count": gp, "G_prime_failure_count": n-gp,
            "G_prime_failure_fraction": (n-gp)/n,
            "G_prime_failure_cp95_lower": gp_l,
            "G_prime_failure_cp95_upper": gp_u,
            "G_prime_only_count": sum(int(r["G_prime"]) and not int(r["G"]) for r in rows),
            "exact_reason_counts": dict(sorted(reason.items())),
            "stop_reason_component_counts": dict(sorted(stops.items())),
            "cyclic_components_total": sum(int(r["cyclic_components"]) for r in rows),
            "largest_complex_max": max(int(r["largest_complex"]) for r in rows),
        })
    else:
        summary["exact_fallback_count"] = sum(int(r["exact_fallback"]) for r in rows)
        summary["solver_status_counts"] = dict(sorted(Counter(r["solver_status"] for r in rows).items()))
    result = {"summary": summary, "integer_moment_sums": integer_sums,
              "cost_interval_width_histogram": dict(sorted(widths.items())),
              "distribution": distribution, "plot_bins": bins, "residual_survival": survival}
    if stage == "bernoulli":
        result["cyclic_components_histogram"] = dict(sorted(Counter(int(r["cyclic_components"]) for r in rows).items()))
        result["largest_complex_histogram"] = dict(sorted(Counter(int(r["largest_complex"]) for r in rows).items()))
        result["cap_affected_draw_counts"] = dict(sorted(Counter(
            key for r in rows for key, count in json.loads(r["stop_summary"]).items() if count > 0
        ).items()))
        result["bound_witness_source_counts"] = dict(sorted(Counter(r["bound_witness_source"] for r in rows).items()))
    return result


def membership_counts(rows, left, right=None):
    definite = possible = 0
    for row in rows:
        low, high = int(row["cost_lower"]), int(row["cost_upper"])
        definite += low >= left and (right is None or high <= right)
        possible += high >= left and (right is None or low <= right)
    return definite, possible


def bernoulli_expectations(rows, common, mixture_info):
    """Implement the frozen analysis protocol from its mathematical definitions."""
    source = common["summary"]
    n, size = source["draws"], source["n"]
    mu = source["poisson_mean"]
    summary = {
        "mixture": source["mixture"], "n": size, "M": n, "mu": mu,
        "poisson_second_moment": source["poisson_raw_second_moment"],
        "g_failures": source["certificate_failure_count"],
        "g_rate": source["certificate_failure_fraction"],
        "g_cp_lo": source["certificate_failure_cp95_lower"],
        "g_cp_hi": source["certificate_failure_cp95_upper"],
        "g_prime_failures": source["G_prime_failure_count"],
        "g_prime_rate": source["G_prime_failure_fraction"],
        "g_prime_cp_lo": source["G_prime_failure_cp95_lower"],
        "g_prime_cp_hi": source["G_prime_failure_cp95_upper"],
        "empirical_cdf_distance_lower": source["empirical_cdf_discrepancy_lower"],
        "empirical_cdf_distance_upper": source["empirical_cdf_discrepancy_upper"],
    }
    unweighted = 2*sum(q*f for q, f in zip(mixture_info["probabilities"], mixture_info["fractions"]))
    summary.update({"lambda_weighted": source["lambda_parameter"], "lambda_unweighted": unweighted,
                    "mu_unweighted": -math.log1p(-unweighted) if unweighted < 1 else None,
                    "exact": source["resolved_count"], "unresolved": source["unresolved_count"],
                    "original_unresolved": sum(r["original_lower"] != r["original_upper"] for r in rows),
                    "newly_resolved": sum(r["original_lower"] != r["original_upper"] and r["cost_lower"] == r["cost_upper"] for r in rows),
                    "original_G_exact": source["G_count"], "weaker_G_prime_exact": source["G_prime_only_count"],
                    "beyond_G_prime_exact": sum(not int(r["G_prime"]) and r["cost_lower"] == r["cost_upper"] for r in rows),
                    "max_residual": source["residual_max"], "max_complex": source["largest_complex_max"],
                    "total_solver_work": source["solver_work_total"],
                    "rescued_count": source["G_prime_only_count"], "rescued_fraction_all": source["G_prime_only_count"]/n,
                    "rescued_fraction_G_failures": source["G_prime_only_count"]/source["G_failure_count"] if source["G_failure_count"] else None})
    for event in ("g", "g_prime"):
        for field in ("rate", "cp_lo", "cp_hi"):
            summary[f"{event}_scaled_{field}"] = size*summary[f"{event}_{field}"]
    for power, prefix in ((1, "mean"), (2, "second_moment")):
        lower = source[f"raw_moment_{power}_lower"]
        upper = source[f"raw_moment_{power}_upper"]
        se_l = source[f"endpoint_raw_moment_{power}_lower_mc_se"]
        se_u = source[f"endpoint_raw_moment_{power}_upper_mc_se"]
        summary.update({prefix+"_lower": lower, prefix+"_upper": upper,
                        prefix+"_endpoint_se_lower": se_l, prefix+"_endpoint_se_upper": se_u,
                        prefix+"_sampling_lo": max(0.0, lower-1.96*se_l),
                        prefix+"_sampling_hi": upper+1.96*se_u,
                        prefix+"_ordinary_se": se_l if source["unresolved_count"] == 0 else None})
    p0_d, p0_p = membership_counts(rows, 0, 0)
    summary.update({"p0_definite_count": p0_d, "p0_possible_count": p0_p,
                    "p0_lower": p0_d/n, "p0_upper": p0_p/n, "poisson_p0": math.exp(-mu),
                    "p0_cp_lo": cp95(p0_d, n)[0], "p0_cp_hi": cp95(p0_p, n)[1]})
    pmf, cdf = poisson_sequence(mu, 100)
    candidates = [k for k in range(1, 101) if all(n*pmf[j] >= 5 for j in range(k)) and n*(1-cdf[k-1]) >= 5]
    cutoff = max(candidates)
    pearson_bins = []
    stat_lower = stat_upper = 0.0
    for k in range(cutoff+1):
        tail = k == cutoff
        definite, possible = membership_counts(rows, k, None if tail else k)
        probability = (1-cdf[k-1]) if tail else pmf[k]
        expectation = n*probability
        closest = min(max(expectation, definite), possible)
        contribution_lower = (closest-expectation)**2/expectation
        contribution_upper = max((definite-expectation)**2, (possible-expectation)**2)/expectation
        stat_lower += contribution_lower
        stat_upper += contribution_upper
        pearson_bins.append({"bin": f"{k}+" if tail else str(k), "definite_count": definite,
                             "possible_count": possible, "expected_count": expectation,
                             "statistic_lower_contribution": contribution_lower,
                             "statistic_upper_contribution": contribution_upper})
    p_l, p_u = chi_square_sf(stat_upper, cutoff), chi_square_sf(stat_lower, cutoff)
    decision = ("reject_for_all_completions" if p_u < 0.05 else
                "fail_to_reject_for_all_completions" if p_l >= 0.05 else "computationally_unresolved")
    summary.update({"chi_square_statistic_lower": stat_lower, "chi_square_statistic_upper": stat_upper,
                    "chi_square_p_lower": p_l, "chi_square_p_upper": p_u,
                    "chi_square_df": cutoff, "chi_square_decision": decision})
    pmf_intervals = []
    for k in range(9):
        tail = k == 8
        definite, possible = membership_counts(rows, k, None if tail else k)
        probability = 1-cdf[7] if tail else pmf[k]
        cp_l, cp_u = cp95(definite, n)[0], cp95(possible, n)[1]
        pmf_intervals.append({"bin": "8+" if tail else str(k), "definite_count": definite,
                              "possible_count": possible, "lower": definite/n, "upper": possible/n,
                              "cp_lo": cp_l, "cp_hi": cp_u, "poisson_probability": probability,
                              "difference_lower": definite/n-probability,
                              "difference_upper": possible/n-probability,
                              "difference_cp_lo": cp_l-probability, "difference_cp_hi": cp_u-probability})
    return {"summary": summary, "pearson_bins": pearson_bins, "pmf_intervals": pmf_intervals,
            "rescued_count": source["G_prime_only_count"],
            "rescued_fraction_all": source["G_prime_only_count"]/n,
            "rescued_fraction_G_failures": source["G_prime_only_count"]/source["G_failure_count"] if source["G_failure_count"] else None,
            "scaled_failure_rates": {key: size*summary[key] for key in
                                     ("g_rate", "g_cp_lo", "g_cp_hi", "g_prime_rate", "g_prime_cp_lo", "g_prime_cp_hi")}}


def local_slopes(expectations, sizes):
    lookup = {(cell["summary"]["mixture"], cell["summary"]["n"]): cell["summary"] for cell in expectations}
    result = []
    spans = list(zip(sizes, sizes[1:])) + [(sizes[0], sizes[2]), (sizes[2], sizes[4])]
    for mixture in sorted({mix for mix, _ in lookup}):
        for event, key in (("g", "g_rate"), ("g_prime", "g_prime_rate")):
            for index, (start, end) in enumerate(spans):
                p_start, p_end = lookup[mixture, start][key], lookup[mixture, end][key]
                slope = math.log(p_end/p_start)/math.log(end/start) if p_start and p_end else None
                result.append({"mixture": mixture, "event": event, "n1": start, "n2": end,
                               "local_log_log_slope": slope,
                               "rate1": p_start, "rate2": p_end,
                               "span": "adjacent" if index < 4 else "prespecified_wide",
                               "status": "defined" if slope is not None else "zero_count_endpoint"})
    return result


def check_records(rows, mixture, size, count, stage, audit):
    label = f"{stage}/{mixture}_n{size}"
    audit.check(len(rows) == count, "cell_row_counts", label+": row count")
    identities = [int(r["replicate"]) for r in rows]
    audit.check(sorted(identities) == list(range(count)), "complete_unique_replicate_identity", label)
    for r in rows:
        tag = label+"/"+r["replicate"]
        lo, hi, residual = (int(r[k]) for k in ("cost_lower", "cost_upper", "residual_size"))
        audit.check(r["mixture"] == mixture and int(r["n"]) == size, "row_cell_identity", tag)
        audit.check(0 <= lo <= hi <= residual <= size, "valid_interval_endpoints", tag)
        pins = [int(value) for value in r["pins"].split(";") if value]
        audit.check(len(pins) == len(set(pins)) == hi and all(0 <= p < size for p in pins),
                    "pin_witness_cardinality_and_range", tag)
        g = int(r["G"] if stage == "bernoulli" else r["certified"])
        audit.check(g in (0, 1) and g == int(int(r["double_rows"]) == 0), "G_double_row_equivalence", tag)
        audit.check(not g or lo == hi, "G_exact_interval", tag)
        if stage == "bernoulli":
            gp = int(r["G_prime"])
            audit.check(gp in (0, 1) and (not g or gp), "G_implies_G_prime", tag)
            audit.check(not gp or lo == hi == int(r["certificate_cycle_count"]), "G_prime_exact_cycle_cost", tag)
            digest = r["network_sha256"]
            audit.check(len(digest) == 64 and all(c in "0123456789abcdef" for c in digest), "network_hash_format", tag)
            stops = json.loads(r["stop_summary"])
            audit.check(all(isinstance(v, int) and v >= 0 for v in stops.values()), "stop_reason_counts", tag)
            audit.check(r["exact_reason"] in {"original_G", "weaker_G_prime", "SCC_exact", "preserved_original_exact",
                                               "intersected_bounds_exact", "bounded_unresolved"}
                        and ((lo < hi) == (r["exact_reason"] == "bounded_unresolved")), "exact_reason_interval_consistency", tag)
            expected_reason = ("original_G" if g else "weaker_G_prime" if gp else
                               "SCC_exact" if r["scc_lower"] == r["scc_upper"] else
                               "preserved_original_exact" if r["original_lower"] == r["original_upper"] else
                               "intersected_bounds_exact" if lo == hi else "bounded_unresolved")
            audit.check(r["exact_reason"] == expected_reason, "exact_reason_derivation", tag)
            audit.check(0 <= int(r["largest_complex"]) <= residual and 0 <= int(r["cyclic_components"]) <= residual
                        and int(r["solver_work"]) >= 0, "SCC_work_summary_ranges", tag)


def compare_baseline_exports(base, cells, audit, root):
    expected = {(c["summary"]["mixture"], c["summary"]["n"]): c for c in cells}
    for filename in ("cell_summary.json", "cell_summary.csv"):
        path = base/filename
        audit.hash(path, root)
        records = json.loads(path.read_text()) if path.suffix == ".json" else read_csv(path)
        keys = [(r["mixture"], int(r["n"])) for r in records]
        audit.check(set(keys) == set(expected) and len(keys) == len(expected), "baseline_export_coverage", str(path))
        for r in records:
            for key, value in expected[(r["mixture"], int(r["n"]))]["summary"].items():
                audit.compare(value, r.get(key), f"{path.name}/{r['mixture']}/{r['n']}/{key}")
            for side in ("lower", "upper"):
                compare_scipy_cp(int(r["certificate_failure_count"]), int(r["draws"]),
                                 r["certificate_failure_cp95_"+side], side, f"baseline/{path.name}/{r['mixture']}/{r['n']}/{side}", audit)
    for filename, section, index in (
        ("cost_distribution.csv", "distribution", "k"),
        ("cost_plot_bins.csv", "plot_bins", "bin"),
        ("residual_survival.csv", "residual_survival", "threshold"),
    ):
        path = base/filename
        audit.hash(path, root)
        rows = read_csv(path)
        lookup = {(mix, size, e[index]): e for (mix, size), cell in expected.items() for e in cell[section]}
        observed = [(r["mixture"], int(r["n"]), int(r[index])) for r in rows]
        audit.check(set(observed) == set(lookup) and len(observed) == len(lookup), "baseline_export_coverage", str(path))
        for row in rows:
            identity = (row["mixture"], int(row["n"]), int(row[index]))
            for key, value in lookup[identity].items():
                audit.compare(value, row[key], f"{path.name}/{identity}/{key}")
            audit.compare(expected[identity[:2]]["summary"]["draws"], row["draws"], f"{path.name}/{identity}/draws")


def compare_table(path, expected_rows, identity_keys, audit, root):
    audit.hash(path, root)
    observed_rows = read_csv(path)
    def identity(row):
        return tuple(str(row[key]) for key in identity_keys)
    lookup = {identity(row): row for row in expected_rows}
    observed_keys = [identity(row) for row in observed_rows]
    audit.check(len(lookup) == len(expected_rows) and len(observed_keys) == len(lookup)
                and set(observed_keys) == set(lookup), "bernoulli_export_coverage", path.name)
    for row in observed_rows:
        key = identity(row)
        if key not in lookup:
            continue
        audit.check(set(row) == set(lookup[key]), "bernoulli_export_field_coverage", f"{path.name}/{key}")
        for field, value in lookup[key].items():
            audit.compare(value, row.get(field), f"{path.name}/{key}/{field}")


def compare_bernoulli_exports(base, cells, expectations, slopes, audit, root, study_protocol_path):
    lookup = {(c["summary"]["mixture"], c["summary"]["n"]): c for c in cells}
    analyses = {(c["summary"]["mixture"], c["summary"]["n"]): c for c in expectations}
    path = base/"summary.json"
    audit.hash(path, root)
    primary = json.loads(path.read_text())
    implementation_path = root/"verification"/"bernoulli_stats.py"
    analysis_path = root/"verification"/"Bernoulli_Analysis_Protocol.json"
    audit.hash(implementation_path, root)
    replay_path = root/"simulation"/"results"/"main"/"run_manifest.json"
    for field, source_path in (("implementation_sha256", implementation_path), ("analysis_protocol_sha256", analysis_path),
                               ("study_protocol_sha256", study_protocol_path), ("replay_manifest_sha256", replay_path)):
        audit.compare(hashlib.sha256(source_path.read_bytes()).hexdigest(), primary[field], f"summary.json/{field}")
    provenance_expected = {f"{mix}_n{size}.csv.gz": common["summary"]["draws"] for (mix, size), common in lookup.items()}
    provenance_observed = [record["file"] for record in primary["provenance"]]
    audit.check(set(provenance_observed) == set(provenance_expected) and len(provenance_observed) == len(provenance_expected),
                "primary_raw_provenance_coverage", "summary.json/provenance")
    for record in primary["provenance"]:
        name = record["file"]
        if name not in provenance_expected:
            continue
        raw_path = replay_path.parent/"raw"/name
        audit.compare(provenance_expected[name], record["rows"], "provenance/"+name+"/rows")
        audit.compare(audit.source_hashes[str(raw_path.relative_to(root))], record["sha256"], "provenance/"+name+"/sha256")
    audit.compare("PASS", primary["status"], "primary statistical status")
    identities = [(c["summary"]["mixture"], c["summary"]["n"]) for c in primary["cells"]]
    audit.check(set(identities) == set(lookup) and len(identities) == len(lookup), "bernoulli_export_coverage", "summary.json/cells")
    for cell in primary["cells"]:
        key = cell["summary"]["mixture"], cell["summary"]["n"]
        if key not in analyses:
            continue
        expected, common = analyses[key], lookup[key]
        audit.check(set(cell["summary"]) == set(expected["summary"]), "bernoulli_export_field_coverage", f"summary.json/{key}")
        for field, value in expected["summary"].items():
            audit.compare(value, cell["summary"].get(field), f"summary.json/{key}/{field}")
        compare_scipy_summary(cell["summary"], f"summary.json/{key}", audit)
        for field in ("exact_reason_counts", "stop_reason_component_counts"):
            audit.compare(common["summary"][field], cell[field], f"summary.json/{key}/{field}")
        pearson = cell["pearson"]
        bins = expected["pearson_bins"]
        audit.compare(len(bins), len(pearson["bins"]), f"summary.json/{key}/Pearson bin count")
        for i, (wanted, actual) in enumerate(zip(bins, pearson["bins"])):
            comparison = {"left": i, "right": None if i == len(bins)-1 else i,
                          "probability": wanted["expected_count"]/expected["summary"]["M"],
                          "expected": wanted["expected_count"], "definite": wanted["definite_count"],
                          "possible": wanted["possible_count"],
                          "contribution_lower": wanted["statistic_lower_contribution"],
                          "contribution_upper": wanted["statistic_upper_contribution"]}
            audit.check(set(comparison) == set(actual), "bernoulli_export_field_coverage", f"Pearson/{key}/{i}")
            for field, value in comparison.items():
                audit.compare(value, actual.get(field), f"Pearson/{key}/{i}/{field}")
        for field in ("df", "statistic_lower", "statistic_upper", "p_lower", "p_upper", "decision"):
            audit.compare(expected["summary"]["chi_square_"+field], pearson[field], f"Pearson/{key}/{field}")
        for side, stat_side in (("lower", "upper"), ("upper", "lower")):
            reference = float(scipy_chi2.sf(pearson["statistic_"+stat_side], pearson["df"]))
            compare_scipy_value(reference, pearson["p_"+side], "SciPy_chi2_sf_export_bound", f"summary.json/{key}/Pearson/{side}", audit)
    total_keys = ("M", "exact", "unresolved", "original_unresolved", "newly_resolved", "original_G_exact",
                  "weaker_G_prime_exact", "beyond_G_prime_exact", "g_failures", "g_prime_failures",
                  "rescued_count", "total_solver_work")
    expected_totals = {key: sum(c["summary"][key] for c in expectations) for key in total_keys}
    for field in ("exact_reason_counts", "stop_reason_component_counts"):
        counts = Counter()
        for cell in cells:
            counts.update(cell["summary"][field])
        expected_totals[field] = dict(sorted(counts.items()))
    expected_totals["pearson_decision_counts"] = dict(sorted(Counter(c["summary"]["chi_square_decision"] for c in expectations).items()))
    audit.check(set(expected_totals) == set(primary["totals"]), "bernoulli_export_field_coverage", "summary.json/totals")
    for field, value in expected_totals.items():
        audit.compare(value, primary["totals"].get(field), f"summary.json/totals/{field}")
    compare_table(base/"cell_summary.csv", [c["summary"] for c in expectations], ("mixture", "n"), audit, root)
    for summary_row in read_csv(base/"cell_summary.csv"):
        compare_scipy_summary(summary_row, f"cell_summary.csv/{summary_row['mixture']}/{summary_row['n']}", audit)
    cdf_rows = []
    pmf_rows = []
    survival_rows = []
    largest_size = max(size for _, size in lookup)
    for key, common in lookup.items():
        mixture, size = key
        n = common["summary"]["draws"]
        for row in common["distribution"]:
            cdf_rows.append({"mixture": mixture, "n": size, **{field: row[field] for field in ("k", "cdf_lower", "cdf_upper", "poisson_cdf")}})
        if size != largest_size:
            continue
        for i, row in enumerate(analyses[key]["pmf_intervals"]):
            pmf_rows.append({"mixture": mixture, "n": size, "M": n, "k": i, "bin_label": row["bin"],
                             "definite_count": row["definite_count"], "possible_count": row["possible_count"],
                             "prob_lower": row["lower"], "prob_upper": row["upper"],
                             "cp_lower": row["cp_lo"], "cp_upper": row["cp_hi"],
                             "poisson_probability": row["poisson_probability"],
                             "difference_lower": row["difference_lower"], "difference_upper": row["difference_upper"],
                             "difference_cp_lower": row["difference_cp_lo"], "difference_cp_upper": row["difference_cp_hi"]})
        for row in common["residual_survival"] + [{"threshold": common["summary"]["residual_max"]+1, "residual_survival_ge": 0.0}]:
            survival_rows.append({"mixture": mixture, "n": size, "M": n, "threshold": row["threshold"],
                                  "surviving_count": round(n*row["residual_survival_ge"]), "survival": row["residual_survival_ge"]})
    compare_table(base/"cdf_intervals.csv", cdf_rows, ("mixture", "n", "k"), audit, root)
    compare_table(base/"pmf_intervals.csv", pmf_rows, ("mixture", "n", "k"), audit, root)
    for row in read_csv(base/"pmf_intervals.csv"):
        for side, count_name in (("lower", "definite_count"), ("upper", "possible_count")):
            compare_scipy_cp(int(row[count_name]), int(row["M"]), row["cp_"+side], side,
                             f"pmf_intervals.csv/{row['mixture']}/{row['k']}/{side}", audit)
    compare_table(base/"residual_survival.csv", survival_rows, ("mixture", "n", "threshold"), audit, root)
    compare_table(base/"local_slopes.csv", slopes, ("mixture", "event", "n1", "n2"), audit, root)
    return {"summary_fields_per_cell": len(expectations[0]["summary"]),
            "cells": len(lookup), "CDF_rows": len(cdf_rows), "PMF_rows": len(pmf_rows),
            "survival_rows": len(survival_rows), "slope_rows": len(slopes), "totals": expected_totals}


def compare_operational(base, cells, expectations, audit, root):
    path = base/"operational_summary.json"
    audit.hash(path, root)
    actual_cells = json.loads(path.read_text())
    lookup = {(c["summary"]["mixture"], c["summary"]["n"]): (c, e) for c, e in zip(cells, expectations)}
    identities = [(c["mixture"], c["n"]) for c in actual_cells]
    audit.check(set(identities) == set(lookup) and len(identities) == len(lookup), "operational_export_coverage", path.name)
    expected_csv_rows = []
    for actual in actual_cells:
        key = actual["mixture"], actual["n"]
        common, statistical = lookup[key]
        a, s = common["summary"], statistical["summary"]
        fields = {"draws": s["M"], "G_failures": s["g_failures"], "G_prime_failures": s["g_prime_failures"],
                  "unresolved": s["unresolved"], "original_unresolved": s["original_unresolved"],
                  "newly_resolved": s["newly_resolved"], "max_residual": s["max_residual"],
                  "max_complex": s["max_complex"], "total_solver_work": s["total_solver_work"],
                  "stop_reason_counts": a["stop_reason_component_counts"], "exact_reason_counts": a["exact_reason_counts"]}
        for field, value in fields.items():
            audit.compare(value, actual.get(field), f"operational/{key}/{field}")
        raw_path = base/"raw"/f"{key[0]}_n{key[1]}.csv.gz"
        audit.compare(audit.source_hashes[str(raw_path.relative_to(root))], actual["raw_sha256"], f"operational/{key}/raw_sha256")
        detail_path = base/"details"/f"{key[0]}_n{key[1]}.jsonl.gz"
        audit.hash(detail_path, root)
        detail_hash = audit.source_hashes[str(detail_path.relative_to(root))]
        audit.compare(detail_hash, actual["detail_sha256"], f"operational/{key}/detail_sha256")
        expected_csv_rows.append({"mixture": key[0], "n": key[1], **fields,
                                  "raw_sha256": audit.source_hashes[str(raw_path.relative_to(root))],
                                  "detail_sha256": detail_hash})
    compare_table(base/"operational_summary.csv", expected_csv_rows, ("mixture", "n"), audit, root)


def compare_numerical_boundaries(base, audit, root):
    path = base/"numerical_boundary_cases.json"
    audit.hash(path, root)
    cases = json.loads(path.read_text())
    audit.compare(hashlib.sha256((root/"verification"/"bernoulli_stats.py").read_bytes()).hexdigest(),
                  cases["implementation_sha256"], "numerical_boundary_cases/implementation_sha256")
    identities = [(r["k"], r["n"]) for r in cases["clopper_pearson"]]
    audit.check(identities == [(k, 10000) for k in range(10001)], "CP_exhaustive_count_coverage", path.name)
    max_cp_error = 0.0
    for row in cases["clopper_pearson"]:
        k, n = row["k"], row["n"]
        reference = [0.0 if k == 0 else float(scipy_beta.ppf(0.025, k, n-k+1)),
                     1.0 if k == n else float(scipy_beta.ppf(0.975, k+1, n-k))]
        for side, wanted in zip(("lower", "upper"), reference):
            error = abs(row[side]-wanted)
            max_cp_error = max(max_cp_error, error)
            compare_scipy_value(wanted, row[side], "SciPy_beta_ppf_all_10001_counts", f"boundary/{k}/{n}/{side}", audit)
    expected_chi = {(df, value) for df in range(1, 21) for value in (0.0, 1e-12, .01, .5, 1.0, 5.0, 10.0, 50.0, 100.0, 1000.0)}
    actual_chi = [(r["df"], r["statistic"]) for r in cases["chi_square"]]
    audit.check(set(actual_chi) == expected_chi and len(actual_chi) == len(expected_chi), "chi_square_numerical_case_coverage", path.name)
    max_chi_error = 0.0
    for row in cases["chi_square"]:
        reference = float(scipy_chi2.sf(row["statistic"], row["df"]))
        max_chi_error = max(max_chi_error, abs(row["survival"]-reference))
        compare_scipy_value(reference, row["survival"], "SciPy_chi2_sf_200_boundary_cases", f"boundary/{row['df']}/{row['statistic']}", audit)
        compare_scipy_value(reference, chi_square_sf(row["statistic"], row["df"]), "independent_gamma_recurrence_vs_SciPy", f"independent/{row['df']}/{row['statistic']}", audit)
    return {"CP_count_cases": len(identities), "CP_endpoints": 2*len(identities),
            "chi_square_cases": len(actual_chi), "CP_maximum_absolute_error": max_cp_error,
            "chi_square_maximum_absolute_error": max_chi_error,
            "relative_tolerance": 2e-10, "absolute_tolerance": 2e-12}


def totals(cells):
    summaries = [cell["summary"] for cell in cells]
    keys = ["draws", "certified_count", "certificate_failure_count", "resolved_count", "unresolved_count",
            "solver_work_total", "G_prime_count", "G_prime_failure_count", "G_prime_only_count", "cyclic_components_total"]
    return {key: sum(s[key] for s in summaries) for key in keys if key in summaries[0]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output or root/"results"/"independent_statistics.json"
    baseline_root = root/"simulation"/"baseline"
    protocol_path = baseline_root/"Protocol_Baseline.json"
    protocol = json.loads(protocol_path.read_text())
    audit = Audit()
    audit.hash(protocol_path, root)
    audit.hash(Path(__file__).resolve(), root)
    expected_cells = [(mix, size) for mix in sorted(protocol["mixtures"]) for size in protocol["sizes"]]
    expected_names = {f"{mix}_n{size}.csv.gz" for mix, size in expected_cells}
    count = protocol["replicates_per_cell"]
    audit.check(len(expected_cells) == protocol["total_cells"], "protocol_internal_consistency", "protocol cell count")
    audit.check(count*len(expected_cells) == protocol["total_networks"], "protocol_internal_consistency", "protocol network count")
    for mixture, info in protocol["mixtures"].items():
        weighted = sum(d*q*f for d, q, f in zip(info["degrees"], info["probabilities"], info["fractions"]))
        audit.compare(info["lambda"], weighted, f"protocol/{mixture}/lambda_formula")
    folders = {"original": baseline_root/"results", "refined": baseline_root/"results"/"refined",
               "bernoulli": root/"simulation"/"results"/"main"}
    have_bernoulli = folders["bernoulli"].joinpath("raw").exists() and {p.name for p in (folders["bernoulli"]/"raw").glob("*.csv.gz")} == expected_names
    if not have_bernoulli:
        parser.error("All 25 Bernoulli raw files are required.")
    stats_base = root/"results"/"bernoulli_statistics"
    stats_names = ("summary.json", "cell_summary.csv", "cdf_intervals.csv", "pmf_intervals.csv", "residual_survival.csv", "local_slopes.csv")
    if have_bernoulli and not all((stats_base/name).is_file() for name in stats_names):
        parser.error("Bernoulli raw files are present but all six primary statistical exports must be ready before independent comparison.")
    if have_bernoulli:
        manifest_path = folders["bernoulli"]/"run_manifest.json"
        audit.hash(manifest_path, root)
        manifest = json.loads(manifest_path.read_text())
        audit.check(manifest["status"] == "completed" and not manifest["pilot"], "Bernoulli_run_complete", str(manifest_path))
        audit.compare(protocol["total_networks"], manifest["total_draws"], "replay manifest total_draws")
        audit.compare(hashlib.sha256(protocol_path.read_bytes()).hexdigest(), manifest["original_protocol_sha256"], "replay manifest protocol hash")
        audit.hash(root/"verification"/"Bernoulli_Analysis_Protocol.json", root)
    stages = ["original", "refined"] + (["bernoulli"] if have_bernoulli else [])
    datasets = {stage: [] for stage in stages}
    identity_digests = {stage: hashlib.sha256() for stage in stages}
    unresolved_expected = {stage: [] for stage in ("original", "refined")}
    transition_counts = Counter()
    expectations = []
    for stage in stages:
        actual_names = {p.name for p in (folders[stage]/"raw").glob("*.csv.gz")}
        audit.check(actual_names == expected_names, "raw_cell_file_coverage", stage)
    for mixture, size in expected_cells:
        stage_rows = {}
        lam = protocol["mixtures"][mixture]["lambda"]
        for stage in stages:
            path = folders[stage]/"raw"/f"{mixture}_n{size}.csv.gz"
            audit.hash(path, root)
            rows = read_csv(path)
            check_records(rows, mixture, size, count, stage, audit)
            rows.sort(key=lambda r: int(r["replicate"]))
            for row in rows:
                identity_digests[stage].update(f"{protocol['master_seed']}:{protocol['mixture_ids'][mixture]}:{size}:{row['replicate']}\n".encode())
            stage_rows[stage] = rows
            datasets[stage].append(aggregate(rows, mixture, size, lam, stage))
            if stage == "bernoulli":
                expectations.append(bernoulli_expectations(rows, datasets[stage][-1], protocol["mixtures"][mixture]))
            if stage in unresolved_expected:
                unresolved_expected[stage].extend(row for row in rows if row["cost_lower"] != row["cost_upper"])
        for original, refined in zip(stage_rows["original"], stage_rows["refined"]):
            label = f"{mixture}/{size}/{original['replicate']}"
            retained = ("mixture", "n", "replicate", "residual_size", "double_rows", "certified", "exact_fallback", "solver_work")
            audit.check(all(original[k] == refined[k] for k in retained), "baseline_refinement_immutable_fields", label)
            olo, ohi, rlo, rhi = (int(row[key]) for row, key in ((original, "cost_lower"), (original, "cost_upper"), (refined, "cost_lower"), (refined, "cost_upper")))
            audit.check(olo <= rlo <= rhi <= ohi, "baseline_refinement_interval_nesting", label)
            if olo == ohi:
                audit.check(original == refined, "baseline_exact_row_preservation", label)
        if have_bernoulli:
            for original, refined, new in zip(stage_rows["original"], stage_rows["refined"], stage_rows["bernoulli"]):
                label = f"{mixture}/{size}/{new['replicate']}"
                audit.check(all(new[k] == original[k] for k in ("mixture", "n", "replicate", "residual_size", "double_rows")), "bernoulli_baseline_immutable_fields", label)
                audit.check(new["G"] == original["certified"], "bernoulli_baseline_G_identity", label)
                lo, hi, rlo, rhi = (int(row[key]) for row, key in ((new,"cost_lower"), (new,"cost_upper"), (refined,"cost_lower"), (refined,"cost_upper")))
                audit.check(int(new["original_lower"]) == rlo and int(new["original_upper"]) == rhi, "bernoulli_recorded_refined_bounds", label)
                slo, shi = int(new["scc_lower"]), int(new["scc_upper"])
                audit.check(0 <= slo <= shi <= int(new["residual_size"]), "bernoulli_SCC_interval_validity", label)
                audit.check(lo == max(slo, rlo) and hi == min(shi, rhi), "bernoulli_interval_intersection_arithmetic", label)
                audit.check(rlo <= lo <= hi <= rhi, "bernoulli_refined_interval_nesting", label)
                if rlo == rhi:
                    audit.check(lo == hi == rlo, "bernoulli_baseline_exact_cost_preservation", label)
                transition_counts[f"baseline_{'exact' if rlo == rhi else 'interval'}_to_bernoulli_{'exact' if lo == hi else 'interval'}"] += 1
    for stage in ("original", "refined"):
        compare_baseline_exports(folders[stage], datasets[stage], audit, root)
        path = folders[stage]/"unresolved_draws.csv"
        audit.hash(path, root)
        audit.check(read_csv(path) == unresolved_expected[stage], "baseline_unresolved_export_rows", stage)
    digests = {stage: value.hexdigest() for stage, value in identity_digests.items()}
    audit.check(len(set(digests.values())) == 1, "seed_identity_digest_equality", "all included stages")
    primary_comparison = None
    numerical_comparison = None
    slopes = []
    if have_bernoulli:
        slopes = local_slopes(expectations, protocol["sizes"])
        primary_comparison = compare_bernoulli_exports(stats_base, datasets["bernoulli"], expectations, slopes, audit, root, protocol_path)
        compare_operational(folders["bernoulli"], datasets["bernoulli"], expectations, audit, root)
        numerical_comparison = compare_numerical_boundaries(stats_base, audit, root)
    cp_validation = verify_cp_identities(audit)
    report = {
        "schema_version": 1,
        "status": "FAIL" if audit.failures else "PASS",
        "scope": "Independent arithmetic aggregation and raw-record integrity; no production modules imported; no network regeneration or dynamic witness verification in this script.",
        "coverage": {"protocol_cells": len(expected_cells), "unique_draws": count*len(expected_cells),
                     "raw_stages": stages, "raw_rows_read": count*len(expected_cells)*len(stages),
                     "Bernoulli_complete": have_bernoulli, "baseline_aggregate_exports_compared": 12,
                     "Bernoulli_statistical_exports_compared": 6 if have_bernoulli else 0},
        "method": {
            "empirical_counts": "Exact integer histograms and sums; all draws retained.",
            "cost_intervals": "Pointwise empirical PMF/CDF and raw moment bounds are partial identification intervals, not confidence intervals. No midpoint imputation.",
            "cdf_discrepancy": "max_k distance from Poisson CDF to [F_lower(k),F_upper(k)] is a conservative lower bound; maximum absolute endpoint distance is an upper bound. All integer k through max upper cost suffice because empirical CDF is then one and the Poisson tail decreases.",
            "endpoint_SE": "sqrt((N sum X^2 - (sum X)^2)/(N^2(N-1))) for X=L^r or U^r; cost moment SE only when every interval is a singleton.",
            "failure_CI": "Pointwise equal-tailed 95% Clopper-Pearson beta quantiles; scipy.special.betaincinv supplies expected limits, checked by independent binomial-tail summation. Every exported primary endpoint is also checked against scipy.stats.beta.ppf.",
            "chi_square_numerics": "Independently implemented integer/half-integer gamma recurrence; every exported p-value bound and 200 numerical cases additionally checked against scipy.stats.chi2.sf.",
            "G_prime_semantics": "Nominal residual-slot graph, retaining duplicate input slots; not the physically essential dependency graph.",
            "seed_identity": "Protocol master_seed, mixture_id, n, and replicate tuples; checks complete identical record identities, not numerical RNG replay.",
            "numeric_comparison_tolerance": {"relative": 2e-10, "absolute": 2e-12},
        },
        "checks": dict(sorted(audit.checks.items())), "failure_count": len(audit.failures),
        "failures": audit.failures, "seed_identity_sha256": digests,
        "transition_counts": dict(sorted(transition_counts.items())),
        "CP_quantile_validation": cp_validation,
        "primary_special_function_validation": numerical_comparison,
        "SciPy_maximum_absolute_errors": dict(sorted(audit.scipy_errors.items())),
        "Bernoulli_primary_export_comparison": primary_comparison,
        "Bernoulli_analysis_expectations": expectations,
        "Bernoulli_local_slopes": slopes,
        "source_sha256": dict(sorted(audit.source_hashes.items())),
        "datasets": {stage: {"totals": totals(cells), "cells": cells} for stage, cells in datasets.items()},
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, sort_keys=True, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"status": report["status"], "coverage": report["coverage"], "failure_count": len(audit.failures),
                      "first_failures": audit.failures[:12], "totals": {stage: totals(cells) for stage, cells in datasets.items()}}, indent=2))
    return 1 if audit.failures else 0


if __name__ == "__main__":
    sys.exit(main())
