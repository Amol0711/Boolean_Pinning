#!/usr/bin/env python3
"""Aggregate every Baseline draw without conditioning on certification or resolution.

Only numpy/scipy and the standard library are required. The default run requires
every predeclared cell and replicate. Partial-identification bounds describe the
observed empirical distribution; they are not confidence intervals for a
population law. Endpoint Monte Carlo standard errors describe only the sampling
variability of the named bound estimators.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import beta, poisson


SCHEMA = {
    "mixture", "n", "replicate", "residual_size", "double_rows", "certified",
    "exact_fallback", "cost_lower", "cost_upper", "solver_status",
    "solver_work", "pins",
}


def boolean(value: str) -> int:
    value = str(value).strip().lower()
    if value in {"1", "true"}:
        return 1
    if value in {"0", "false"}:
        return 0
    raise ValueError(f"Invalid Boolean value {value!r}")


def cp_interval(count: int, total: int, alpha: float = 0.05) -> tuple[float, float]:
    """Two-sided equal-tailed Clopper--Pearson interval, including endpoints."""
    if not 0 <= count <= total or total < 1:
        raise ValueError("Invalid binomial count or denominator")
    lo = 0.0 if count == 0 else float(beta.ppf(alpha / 2, count, total - count + 1))
    hi = 1.0 if count == total else float(beta.ppf(1 - alpha / 2, count + 1, total - count))
    return lo, hi


def mean_se(values: np.ndarray) -> float | None:
    return float(np.std(values, ddof=1) / math.sqrt(len(values))) if len(values) > 1 else None


def empirical_bounds(lower: np.ndarray, upper: np.ndarray, mu: float) -> dict:
    """Pointwise bounds and valid conservative bounds for empirical Kolmogorov distance."""
    if len(lower) != len(upper) or len(lower) == 0 or np.any(lower < 0) or np.any(lower > upper):
        raise ValueError("Invalid integer cost intervals")
    total = len(lower)
    maximum = int(upper.max(initial=0))
    grid = np.arange(maximum + 1, dtype=np.int64)
    low_hist = np.bincount(lower, minlength=maximum + 1)
    high_hist = np.bincount(upper, minlength=maximum + 1)
    cdf_lo = np.cumsum(high_hist) / total
    cdf_hi = np.cumsum(low_hist) / total
    resolved = lower == upper
    pmf_lo = np.bincount(lower[resolved], minlength=maximum + 1) / total
    # Intervals starting at k count from k onward; intervals ending at k cease
    # to count at k+1. This also includes every unresolved draw.
    pmf_hi = (np.cumsum(low_hist) - np.r_[0, np.cumsum(high_hist)[:-1]]) / total
    target_cdf = poisson.cdf(grid, mu)
    distance_lo = float(np.maximum.reduce((cdf_lo - target_cdf, target_cdf - cdf_hi,
                                          np.zeros_like(target_cdf))).max(initial=0))
    distance_hi = float(np.maximum(np.abs(cdf_lo - target_cdf),
                                   np.abs(cdf_hi - target_cdf)).max(initial=0))
    # Below zero the discrepancy is zero. Above maximum all possible empirical
    # CDFs are one, so their discrepancy decreases monotonically. Hence this
    # finite grid computes bounds for the supremum over ALL integer thresholds.
    return dict(grid=grid, cdf_lo=cdf_lo, cdf_hi=cdf_hi, pmf_lo=pmf_lo,
                pmf_hi=pmf_hi, target_cdf=target_cdf,
                target_pmf=poisson.pmf(grid, mu),
                distance_lo=distance_lo, distance_hi=distance_hi)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Use an explicit schema for empty output {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_cells(raw_dir: Path, protocol: dict) -> tuple[dict, list[dict], list[dict]]:
    expected = {(m, int(n)) for m in protocol["mixtures"] for n in protocol["sizes"]}
    per_cell = int(protocol["replicates_per_cell"])
    cells = defaultdict(list)
    fingerprints = []
    unresolved = []
    paths = sorted(raw_dir.glob("*.csv.gz"))
    if not paths:
        raise FileNotFoundError(f"No raw CSV.gz files in {raw_dir}")
    for path in paths:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        fingerprints.append({"file": path.name, "sha256": digest})
        with gzip.open(path, "rt", newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if not SCHEMA.issubset(set(reader.fieldnames or [])):
                raise ValueError(f"Missing columns in {path}: {SCHEMA - set(reader.fieldnames or [])}")
            for line, row in enumerate(reader, 2):
                try:
                    m = row["mixture"]
                    values = {key: int(row[key]) for key in (
                        "n", "replicate", "residual_size", "double_rows", "cost_lower",
                        "cost_upper", "solver_work")}
                    values.update(mixture=m, certified=boolean(row["certified"]),
                                  exact_fallback=boolean(row["exact_fallback"]),
                                  solver_status=row["solver_status"], pins=row["pins"])
                    key = (m, values["n"])
                    if key not in expected:
                        raise ValueError(f"Undeclared cell {key}")
                    r = values["residual_size"]
                    lo, hi = values["cost_lower"], values["cost_upper"]
                    if not 0 <= values["replicate"] < per_cell:
                        raise ValueError("Replicate index out of protocol range")
                    if not 0 <= lo <= hi <= r <= values["n"]:
                        raise ValueError("Invalid residual size or cost interval")
                    if not 0 <= values["double_rows"] <= r or values["solver_work"] < 0:
                        raise ValueError("Invalid double-row or solver-work count")
                    if values["certified"] != (values["double_rows"] == 0):
                        raise ValueError("Certificate must equal the zero-double-row event")
                    if values["certified"] and lo != hi:
                        raise ValueError("Certified cost is not resolved")
                    if values["certified"] and values["exact_fallback"]:
                        raise ValueError("Certified draw should not use the fallback solver")
                    if values["exact_fallback"] and lo != hi:
                        raise ValueError("exact_fallback=1 denotes successful exact resolution")
                    if not values["solver_status"]:
                        raise ValueError("Empty solver status")
                    pins = [] if not values["pins"] else [int(s) for s in values["pins"].split(";")]
                    if len(set(pins)) != len(pins) or any(not 0 <= p < values["n"] for p in pins):
                        raise ValueError("Invalid returned pin set")
                    if lo == hi and len(pins) != hi:
                        raise ValueError("Resolved optimum and returned pin count differ")
                    cells[key].append(values)
                    if lo != hi:
                        unresolved.append(values.copy())
                except (ValueError, KeyError) as error:
                    raise ValueError(f"{path.name}:{line}: {error}") from error
    if set(cells) != expected:
        raise ValueError(f"Missing cells: {sorted(expected - set(cells))}")
    for key, rows in cells.items():
        ids = [row["replicate"] for row in rows]
        if len(ids) != per_cell or len(set(ids)) != per_cell:
            raise ValueError(f"Incomplete or duplicated replicates in {key}: {len(ids)} rows, {len(set(ids))} unique")
        rows.sort(key=lambda row: row["replicate"])
    total = sum(map(len, cells.values()))
    if total != int(protocol["total_networks"]) or len(cells) != int(protocol["total_cells"]):
        raise ValueError("Protocol totals do not match complete raw records")
    return cells, fingerprints, unresolved


def summarize(cells: dict, protocol: dict) -> tuple[list, list, list, list]:
    summary, distributions, survival, plot_bins = [], [], [], []
    for (mixture, n), rows in sorted(cells.items()):
        total = len(rows)
        lower = np.array([x["cost_lower"] for x in rows], dtype=np.int64)
        upper = np.array([x["cost_upper"] for x in rows], dtype=np.int64)
        residual = np.array([x["residual_size"] for x in rows], dtype=np.int64)
        cert_count = sum(x["certified"] for x in rows)
        failures = total - cert_count
        unresolved_count = int(np.count_nonzero(lower != upper))
        resolved_count = total - unresolved_count
        fallback_count = sum(x["exact_fallback"] for x in rows)
        if cert_count + fallback_count != resolved_count:
            raise ValueError(f"Resolved accounting inconsistent in {(mixture,n)}")
        lam = float(protocol["mixtures"][mixture]["lambda"])
        mu = -math.log1p(-lam)
        cp_lo, cp_hi = cp_interval(failures, total)
        bounds = empirical_bounds(lower, upper, mu)
        item = dict(mixture=mixture, n=n, draws=total, lambda_parameter=lam,
                    poisson_mean=mu, poisson_raw_second_moment=mu + mu*mu,
                    certified_count=cert_count, certificate_failure_count=failures,
                    certificate_failure_fraction=failures/total,
                    certificate_failure_cp95_lower=cp_lo,
                    certificate_failure_cp95_upper=cp_hi,
                    exact_fallback_count=fallback_count,
                    resolved_count=resolved_count, unresolved_count=unresolved_count,
                    resolved_fraction=resolved_count/total,
                    residual_mean=float(np.mean(residual)),
                    residual_mean_mc_se=mean_se(residual.astype(float)),
                    residual_max=int(residual.max()),
                    double_rows_total=sum(x["double_rows"] for x in rows),
                    solver_work_total=sum(x["solver_work"] for x in rows),
                    solver_work_max=max(x["solver_work"] for x in rows),
                    empirical_cdf_discrepancy_lower=bounds["distance_lo"],
                    empirical_cdf_discrepancy_upper=bounds["distance_hi"])
        for power in (1, 2):
            lo_moment = lower.astype(float) ** power
            hi_moment = upper.astype(float) ** power
            item[f"raw_moment_{power}_lower"] = float(np.mean(lo_moment))
            item[f"raw_moment_{power}_upper"] = float(np.mean(hi_moment))
            item[f"endpoint_raw_moment_{power}_lower_mc_se"] = mean_se(lo_moment)
            item[f"endpoint_raw_moment_{power}_upper_mc_se"] = mean_se(hi_moment)
            item[f"resolved_raw_moment_{power}_mc_se"] = mean_se(lo_moment) if not unresolved_count else None
        item["solver_status_counts"] = json.dumps(dict(sorted(Counter(x["solver_status"] for x in rows).items())), sort_keys=True)
        summary.append(item)
        for i, k in enumerate(bounds["grid"]):
            distributions.append(dict(mixture=mixture, n=n, draws=total, k=int(k),
                                      pmf_lower=float(bounds["pmf_lo"][i]), pmf_upper=float(bounds["pmf_hi"][i]),
                                      cdf_lower=float(bounds["cdf_lo"][i]), cdf_upper=float(bounds["cdf_hi"][i]),
                                      poisson_pmf=float(bounds["target_pmf"][i]), poisson_cdf=float(bounds["target_cdf"][i])))
        # k=11 is a single >10 tail bin. The denominator remains all draws.
        for k in range(12):
            if k < 11:
                certainly = (lower == k) & (upper == k)
                possibly = (lower <= k) & (upper >= k)
                target = float(poisson.pmf(k, mu))
            else:
                certainly, possibly = lower > 10, upper > 10
                target = float(poisson.sf(10, mu))
            plot_bins.append(dict(mixture=mixture, n=n, bin=k, draws=total,
                                  pmf_lower=float(np.mean(certainly)), pmf_upper=float(np.mean(possibly)),
                                  poisson_pmf=target))
        histogram = np.bincount(residual)
        survivals = np.cumsum(histogram[::-1])[::-1] / total
        for r, probability in enumerate(survivals):
            survival.append(dict(mixture=mixture, n=n, draws=total, threshold=r,
                                 residual_survival_ge=float(probability)))
    return summary, distributions, survival, plot_bins


def self_test() -> None:
    from itertools import product
    from tempfile import TemporaryDirectory
    lo = np.array([0, 1, 0], dtype=np.int64)
    hi = np.array([0, 1, 2], dtype=np.int64)
    bounds = empirical_bounds(lo, hi, 0.7)
    np.testing.assert_allclose(bounds["pmf_lo"], [1/3, 1/3, 0])
    np.testing.assert_allclose(bounds["pmf_hi"], [2/3, 2/3, 1/3])
    np.testing.assert_allclose(bounds["cdf_lo"], [1/3, 2/3, 1])
    np.testing.assert_allclose(bounds["cdf_hi"], [2/3, 1, 1])
    for completion in product(*(range(int(a), int(b)+1) for a, b in zip(lo, hi))):
        exact = empirical_bounds(np.array(completion), np.array(completion), 0.7)
        assert bounds["distance_lo"] - 1e-14 <= exact["distance_lo"] <= bounds["distance_hi"] + 1e-14
        assert math.isclose(exact["distance_lo"], exact["distance_hi"])
    cp0 = cp_interval(0, 10000)
    assert cp0[0] == 0 and math.isclose(cp0[1], 1 - 0.025**(1/10000), abs_tol=1e-14)
    cpall = cp_interval(10000, 10000)
    assert cpall[1] == 1 and math.isclose(cpall[0], 0.025**(1/10000), abs_tol=1e-14)
    empty_residual = empirical_bounds(np.array([0, 0]), np.array([0, 0]), 1.0)
    assert math.isclose(empty_residual["distance_lo"], 1-math.exp(-1))
    # Exercise ingestion/accounting with an explicitly synthetic three-draw
    # fixture: one certified draw, one exact fallback, and one unresolved draw.
    toy_protocol = {"mixtures": {"T": {"lambda": 0.4}}, "sizes": [3],
                    "replicates_per_cell": 3, "total_cells": 1, "total_networks": 3}
    toy_rows = [
        dict(mixture="T", n=3, replicate=0, residual_size=0, double_rows=0,
             certified=1, exact_fallback=0, cost_lower=0, cost_upper=0,
             solver_status="certified", solver_work=0, pins=""),
        dict(mixture="T", n=3, replicate=1, residual_size=2, double_rows=1,
             certified=0, exact_fallback=1, cost_lower=1, cost_upper=1,
             solver_status="exact", solver_work=10, pins="1"),
        dict(mixture="T", n=3, replicate=2, residual_size=2, double_rows=1,
             certified=0, exact_fallback=0, cost_lower=0, cost_upper=2,
             solver_status="work_cap", solver_work=10, pins="0;1")]
    with TemporaryDirectory(prefix="boolean_baseline_statistics_test_") as tmp:
        fixture = Path(tmp) / "synthetic.csv.gz"
        def save_fixture():
            with gzip.open(fixture, "wt", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(toy_rows[0]))
                writer.writeheader()
                writer.writerows(toy_rows)
        save_fixture()
        loaded, _, unresolved = load_cells(Path(tmp), toy_protocol)
        summary, _, _, _ = summarize(loaded, toy_protocol)
        assert len(unresolved) == 1 and summary[0]["draws"] == 3
        assert summary[0]["certified_count"] == summary[0]["exact_fallback_count"] == summary[0]["unresolved_count"] == 1
        assert summary[0]["resolved_raw_moment_1_mc_se"] is None
        assert math.isclose(summary[0]["raw_moment_1_lower"], 1/3)
        assert math.isclose(summary[0]["raw_moment_1_upper"], 1)
        toy_rows[-1]["replicate"] = 1
        save_fixture()
        try:
            load_cells(Path(tmp), toy_protocol)
        except ValueError as error:
            assert "duplicated" in str(error)
        else:
            raise AssertionError("Duplicate/missing replicate fixture was accepted")
    print("PASS: interval bounds, exhaustive completions, complete-resolution identity, CP endpoints, zero-cost tail, all-draw accounting, unresolved retention, and duplicate/missing-replicate rejection.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--protocol", type=Path)
    parser.add_argument("--raw-dir", type=Path)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    protocol_path = args.protocol or args.root / "Protocol_Baseline.json"
    protocol = json.loads(protocol_path.read_text())
    raw_dir = args.raw_dir or args.root / "results" / "raw"
    out_dir = args.out_dir or args.root / "results"
    cells, fingerprints, unresolved = load_cells(raw_dir, protocol)
    summary, distributions, survival, plot_bins = summarize(cells, protocol)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(out_dir / "cell_summary.csv", summary)
    (out_dir / "cell_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    write_csv(out_dir / "cost_distribution.csv", distributions)
    write_csv(out_dir / "cost_plot_bins.csv", plot_bins)
    write_csv(out_dir / "residual_survival.csv", survival)
    with (out_dir / "unresolved_draws.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["mixture", "n", "replicate", "residual_size", "double_rows",
                                                   "certified", "exact_fallback", "cost_lower", "cost_upper",
                                                   "solver_status", "solver_work", "pins"])
        writer.writeheader()
        writer.writerows(sorted(unresolved, key=lambda x: (x["mixture"], x["n"], x["replicate"])))
    audit = dict(
        status="PASS", all_protocol_cells_present=True, all_replicates_present_once=True,
        total_cells=len(cells), total_draws=sum(x["draws"] for x in summary),
        total_certified=sum(x["certified_count"] for x in summary),
        total_uncertified=sum(x["certificate_failure_count"] for x in summary),
        total_exact_fallback=sum(x["exact_fallback_count"] for x in summary),
        total_unresolved=len(unresolved), exclusions=0,
        protocol_sha256=hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
        raw_file_sha256=fingerprints,
        bounds_interpretation="All-draw partial-identification bounds for the observed empirical law, not population confidence intervals.",
        cdf_discrepancy_interpretation="Valid conservative bounds for the supremum discrepancy between any compatible completed empirical CDF and the specified Poisson CDF over all integer thresholds. These are not population-distance confidence intervals.",
        mc_se_interpretation="Endpoint SE columns are standard errors of the named lower/upper sample-moment estimators only. A cost-moment SE is populated only when every draw in the cell is resolved.",
        cp_interpretation="Pointwise two-sided 95% Clopper-Pearson intervals for each certificate-failure probability; no simultaneous coverage claim.")
    (out_dir / "statistical_audit.json").write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: audit[key] for key in ("status", "total_cells", "total_draws", "total_uncertified", "total_unresolved", "exclusions")}, indent=2))


if __name__ == "__main__":
    main()
