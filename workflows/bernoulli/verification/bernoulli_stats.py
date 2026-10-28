#!/usr/bin/env python3
"""Complete-sample Bernoulli statistical summaries; standard library, no solver imports.

Methods are fixed in Bernoulli_Analysis_Protocol.json. Computational intervals and
sampling uncertainty are separate outputs. No unresolved cost is imputed.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ALPHA = 0.05
Z95 = 1.96


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def beta_fraction(a, b, x):
    """Continued fraction in the regularized incomplete beta function."""
    floor = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < floor: d = floor
    d = 1.0 / d
    result = d
    for m in range(1, 10001):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < floor: d = floor
        c = 1.0 + aa / c
        if abs(c) < floor: c = floor
        d = 1.0 / d
        result *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < floor: d = floor
        c = 1.0 + aa / c
        if abs(c) < floor: c = floor
        d = 1.0 / d
        change = d * c
        result *= change
        if abs(change - 1.0) <= 3e-14:
            return result
    raise ArithmeticError("Incomplete beta fraction did not converge")


def regularized_beta(x, a, b):
    if x <= 0.0: return 0.0
    if x >= 1.0: return 1.0
    factor = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                      + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        result = factor * beta_fraction(a, b, x) / a
    else:
        result = 1.0 - factor * beta_fraction(b, a, 1.0 - x) / b
    return min(1.0, max(0.0, result))


def beta_quantile(p, a, b):
    lo, hi = 0.0, 1.0
    for _ in range(65):
        mid = (lo + hi) / 2
        if regularized_beta(mid, a, b) < p: lo = mid
        else: hi = mid
    return (lo + hi) / 2


def clopper_pearson(k, n, alpha=ALPHA):
    assert 0 <= k <= n and n > 0
    lo = 0.0 if k == 0 else beta_quantile(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta_quantile(1 - alpha / 2, k + 1, n - k)
    return lo, hi


def chi_square_sf(x, df):
    """Chi-square survival for integer degrees of freedom via gamma recurrence."""
    assert df >= 1 and isinstance(df, int) and x >= 0
    if x == 0: return 1.0
    if math.isinf(x): return 0.0
    y = x / 2
    if df % 2 == 0:
        shape, result = 1.0, math.exp(-y)
    else:
        shape, result = 0.5, math.erfc(math.sqrt(y))
    while shape < df / 2:
        result += math.exp(shape * math.log(y) - y - math.lgamma(shape + 1))
        shape += 1
    return min(1.0, max(0.0, result))


def moment_estimates(values):
    n = len(values)
    assert n > 1
    total = sum(values)
    squares = sum(x * x for x in values)
    # Integer arithmetic prevents cancellation before conversion to floating point.
    se = math.sqrt((n * squares - total * total) / (n * n * (n - 1)))
    return total / n, se


def poisson_probabilities(mu, upto):
    result = [math.exp(-mu)]
    for k in range(1, upto + 1): result.append(result[-1] * mu / k)
    return result


def membership_counts(lower, upper, left, right=None):
    definite = sum(l >= left and (right is None or u <= right)
                   for l, u in zip(lower, upper))
    possible = sum(u >= left and (right is None or l <= right)
                   for l, u in zip(lower, upper))
    return definite, possible


def poisson_bins(mu, n):
    pmf = poisson_probabilities(mu, 100)
    candidates = []
    for k in range(1, len(pmf)):
        probabilities = pmf[:k] + [max(0.0, 1.0 - math.fsum(pmf[:k]))]
        if min(probabilities) * n >= 5.0: candidates.append((k, probabilities))
    if not candidates: raise ValueError("No admissible fixed-null Pearson bins")
    k, probabilities = candidates[-1]
    return [(j, j, probabilities[j]) for j in range(k)] + [(k, None, probabilities[-1])]


def pearson_bounds(lower, upper, mu):
    n = len(lower)
    records, stat_lower, stat_upper = [], 0.0, 0.0
    for left, right, probability in poisson_bins(mu, n):
        definite, possible = membership_counts(lower, upper, left, right)
        expected = n * probability
        gap = max(definite - expected, expected - possible, 0.0)
        term_lo = gap * gap / expected
        term_hi = max((definite - expected) ** 2, (possible - expected) ** 2) / expected
        stat_lower += term_lo; stat_upper += term_hi
        records.append(dict(left=left, right=right, probability=probability,
                            expected=expected, definite=definite, possible=possible,
                            contribution_lower=term_lo, contribution_upper=term_hi))
    df = len(records) - 1
    p_lo, p_hi = chi_square_sf(stat_upper, df), chi_square_sf(stat_lower, df)
    decision = ('reject_for_all_completions' if p_hi < ALPHA else
                'fail_to_reject_for_all_completions' if p_lo >= ALPHA else
                'computationally_unresolved')
    return dict(bins=records, df=df, statistic_lower=stat_lower,
                statistic_upper=stat_upper, p_lower=p_lo, p_upper=p_hi, decision=decision)


def summarize_cell(records, mixture_info):
    mixture, n = records[0]['mixture'], int(records[0]['n'])
    m = len(records)
    lower = [int(r['cost_lower']) for r in records]
    upper = [int(r['cost_upper']) for r in records]
    oldlo = [int(r['original_lower']) for r in records]
    oldhi = [int(r['original_upper']) for r in records]
    g = [int(r['G']) for r in records]
    gp = [int(r['G_prime']) for r in records]
    for i, r in enumerate(records):
        assert 0 <= oldlo[i] <= lower[i] <= upper[i] <= oldhi[i] <= int(r['residual_size'])
        assert g[i] in (0, 1) and gp[i] in (0, 1) and g[i] <= gp[i]
        if gp[i]: assert lower[i] == upper[i] == int(r['certificate_cycle_count'])
        pins = [] if not r['pins'] else r['pins'].split(';')
        assert len(pins) == len(set(pins)) == upper[i]
    lam = mixture_info['lambda']
    mu = -math.log1p(-lam)
    unweighted = 2 * sum(q * f for q, f in zip(mixture_info['probabilities'], mixture_info['fractions']))
    s = dict(mixture=mixture, n=n, M=m, lambda_weighted=lam, mu=mu,
             lambda_unweighted=unweighted,
             mu_unweighted=-math.log1p(-unweighted) if unweighted < 1 else None,
             poisson_second_moment=mu + mu * mu,
             exact=sum(l == u for l, u in zip(lower, upper)),
             unresolved=sum(l != u for l, u in zip(lower, upper)),
             original_unresolved=sum(l != u for l, u in zip(oldlo, oldhi)),
             newly_resolved=sum(ol != ou and l == u for ol, ou, l, u in zip(oldlo, oldhi, lower, upper)),
             original_G_exact=sum(g), weaker_G_prime_exact=sum(gp) - sum(g),
             beyond_G_prime_exact=sum(not flag and l == u for flag, l, u in zip(gp, lower, upper)),
             max_residual=max(int(r['residual_size']) for r in records),
             max_complex=max(int(r['largest_complex']) for r in records),
             total_solver_work=sum(int(r['solver_work']) for r in records))
    for key, flags in [('g', g), ('g_prime', gp)]:
        failures = m - sum(flags)
        lo, hi = clopper_pearson(failures, m)
        s.update({key + '_failures': failures, key + '_rate': failures / m,
                  key + '_cp_lo': lo, key + '_cp_hi': hi,
                  key + '_scaled_rate': n * failures / m,
                  key + '_scaled_cp_lo': n * lo, key + '_scaled_cp_hi': n * hi})
    s['rescued_count'] = sum(gp) - sum(g)
    s['rescued_fraction_all'] = s['rescued_count'] / m
    s['rescued_fraction_G_failures'] = (s['rescued_count'] / s['g_failures']
                                       if s['g_failures'] else None)
    for power, label in [(1, 'mean'), (2, 'second_moment')]:
        ml, sl = moment_estimates([v ** power for v in lower])
        mu_, su = moment_estimates([v ** power for v in upper])
        s.update({label + '_lower': ml, label + '_upper': mu_,
                  label + '_endpoint_se_lower': sl, label + '_endpoint_se_upper': su,
                  label + '_sampling_lo': max(0.0, ml - Z95 * sl),
                  label + '_sampling_hi': mu_ + Z95 * su,
                  label + '_ordinary_se': sl if s['unresolved'] == 0 else None})
    definite, possible = membership_counts(lower, upper, 0, 0)
    s.update(p0_definite_count=definite, p0_possible_count=possible,
             p0_lower=definite / m, p0_upper=possible / m,
             p0_cp_lo=clopper_pearson(definite, m)[0],
             p0_cp_hi=clopper_pearson(possible, m)[1], poisson_p0=math.exp(-mu))
    pmf = poisson_probabilities(mu, max(upper))
    cdf_records, dl, du = [], 0.0, 0.0
    cdf_poisson = 0.0
    for k, probability in enumerate(pmf):
        cdf_poisson = min(1.0, cdf_poisson + probability)
        fl, fu = sum(v <= k for v in upper) / m, sum(v <= k for v in lower) / m
        dl = max(dl, fl - cdf_poisson, cdf_poisson - fu)
        du = max(du, abs(fl - cdf_poisson), abs(fu - cdf_poisson))
        cdf_records.append(dict(mixture=mixture, n=n, k=k, cdf_lower=fl,
                                cdf_upper=fu, poisson_cdf=cdf_poisson))
    s.update(empirical_cdf_distance_lower=dl, empirical_cdf_distance_upper=du)
    pearson = pearson_bounds(lower, upper, mu)
    for key in ['df', 'statistic_lower', 'statistic_upper', 'p_lower', 'p_upper', 'decision']:
        s['chi_square_' + key] = pearson[key]
    reasons = Counter(r['exact_reason'] for r in records)
    stops = Counter()
    for r in records: stops.update(json.loads(r['stop_summary']))
    detail = dict(summary=s, exact_reason_counts=dict(sorted(reasons.items())),
                  stop_reason_component_counts=dict(sorted(stops.items())), pearson=pearson)
    pmf_records, survival = [], []
    if n == 25600:
        pmf8 = poisson_probabilities(mu, 7)
        for k in range(9):
            right = None if k == 8 else k
            definite, possible = membership_counts(lower, upper, k, right)
            pl, pu = definite / m, possible / m
            cp_lo, cp_hi = clopper_pearson(definite, m)[0], clopper_pearson(possible, m)[1]
            theoretical = (max(0.0, 1.0 - math.fsum(pmf8)) if k == 8 else pmf8[k])
            pmf_records.append(dict(mixture=mixture, n=n, M=m, k=k,
                                    bin_label='8+' if k == 8 else str(k),
                                    definite_count=definite, possible_count=possible,
                                    prob_lower=pl, prob_upper=pu, cp_lower=cp_lo,
                                    cp_upper=cp_hi, poisson_probability=theoretical,
                                    difference_lower=pl-theoretical, difference_upper=pu-theoretical,
                                    difference_cp_lower=cp_lo-theoretical, difference_cp_upper=cp_hi-theoretical))
        residuals = [int(r['residual_size']) for r in records]
        for k in range(max(residuals) + 2):
            count = sum(v >= k for v in residuals)
            survival.append(dict(mixture=mixture, n=n, M=m, threshold=k,
                                 surviving_count=count, survival=count / m))
    return detail, cdf_records, pmf_records, survival


def write_csv(path, rows):
    assert rows
    with Path(path).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def special_function_cases(output):
    """Export the entire 10,000-trial CP count range for an external oracle."""
    records = []
    for k in range(10001):
        lo, hi = clopper_pearson(k, 10000)
        records.append(dict(k=k, n=10000, lower=lo, upper=hi))
    chi = [dict(df=df, statistic=x, survival=chi_square_sf(x, df))
           for df in range(1, 21)
           for x in [0.0, 1e-12, 0.01, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 1000.0]]
    output.mkdir(parents=True, exist_ok=True)
    result = dict(schema_version='1.0', clopper_pearson=records, chi_square=chi,
                  implementation_sha256=digest(__file__))
    (output/'numerical_boundary_cases.json').write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print(json.dumps(dict(cp_cases=len(records), chi_square_cases=len(chi))))


def run(raw_dir, output, study_protocol):
    protocol = json.loads(study_protocol.read_text())
    manifest_path = raw_dir.parent / 'run_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    assert manifest['status'] == 'completed' and not manifest['pilot']
    assert manifest['total_draws'] == 250000
    assert manifest['original_protocol_sha256'] == digest(study_protocol)
    m = protocol['replicates_per_cell']
    expected = [(mixture, n) for mixture in protocol['mixtures'] for n in protocol['sizes']]
    assert len(expected) == 25 and m == 10000
    cells, cdfs, pmfs, survivals, provenance = [], [], [], [], []
    for mixture, n in expected:
        path = raw_dir / f'{mixture}_n{n}.csv.gz'
        with gzip.open(path, 'rt', newline='') as f: records = list(csv.DictReader(f))
        assert len(records) == m
        assert {int(r['replicate']) for r in records} == set(range(m))
        assert all(r['mixture'] == mixture and int(r['n']) == n for r in records)
        records.sort(key=lambda r: int(r['replicate']))
        cell, cdf, pmf, survival = summarize_cell(records, protocol['mixtures'][mixture])
        cells.append(cell); cdfs.extend(cdf); pmfs.extend(pmf); survivals.extend(survival)
        provenance.append(dict(file=path.name, sha256=digest(path), rows=m))
    slopes = []
    index = {(c['summary']['mixture'], c['summary']['n']): c['summary'] for c in cells}
    pairs = list(zip(protocol['sizes'], protocol['sizes'][1:])) + [(100, 1600), (1600, 25600)]
    for mixture in protocol['mixtures']:
        for pair_index, (n1, n2) in enumerate(pairs):
            for event in ['g', 'g_prime']:
                p1, p2 = index[mixture, n1][event+'_rate'], index[mixture, n2][event+'_rate']
                slope = math.log(p2/p1) / math.log(n2/n1) if p1 and p2 else None
                slopes.append(dict(mixture=mixture, event=event, n1=n1, n2=n2,
                                   span='adjacent' if pair_index < 4 else 'prespecified_wide',
                                   rate1=p1, rate2=p2, local_log_log_slope=slope,
                                   status='defined' if slope is not None else 'zero_count_endpoint'))
    totals = {}
    for key in ['M', 'exact', 'unresolved', 'original_unresolved', 'newly_resolved',
                'original_G_exact', 'weaker_G_prime_exact', 'beyond_G_prime_exact',
                'g_failures', 'g_prime_failures', 'rescued_count', 'total_solver_work']:
        totals[key] = sum(c['summary'][key] for c in cells)
    assert totals['M'] == 250000
    assert totals['original_G_exact'] + totals['weaker_G_prime_exact'] + totals['beyond_G_prime_exact'] + totals['unresolved'] == 250000
    exact_reasons, stop_reasons = Counter(), Counter()
    for c in cells:
        exact_reasons.update(c['exact_reason_counts']); stop_reasons.update(c['stop_reason_component_counts'])
    totals['exact_reason_counts'] = dict(sorted(exact_reasons.items()))
    totals['stop_reason_component_counts'] = dict(sorted(stop_reasons.items()))
    totals['pearson_decision_counts'] = dict(sorted(Counter(c['summary']['chi_square_decision'] for c in cells).items()))
    output.mkdir(parents=True, exist_ok=True)
    summary = dict(schema_version='1.0', status='PASS', cells=cells, totals=totals,
                   provenance=provenance, study_protocol_sha256=digest(study_protocol),
                   replay_manifest_sha256=digest(manifest_path),
                   analysis_protocol_sha256=digest(Path(__file__).with_name('Bernoulli_Analysis_Protocol.json')),
                   implementation_sha256=digest(__file__))
    (output/'summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False)+'\n')
    write_csv(output/'cell_summary.csv', [c['summary'] for c in cells])
    write_csv(output/'cdf_intervals.csv', cdfs)
    write_csv(output/'pmf_intervals.csv', pmfs)
    write_csv(output/'residual_survival.csv', survivals)
    write_csv(output/'local_slopes.csv', slopes)
    print(json.dumps(dict(status='PASS', cells=len(cells), totals=totals), sort_keys=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir', type=Path, default=ROOT/'simulation/results/main/raw')
    parser.add_argument('--output', type=Path, default=ROOT/'results/bernoulli_statistics')
    parser.add_argument('--study-protocol', type=Path, default=ROOT/'simulation/baseline/Protocol_Baseline.json')
    parser.add_argument('--special-function-checks-only', action='store_true')
    args = parser.parse_args()
    if args.special_function_checks_only: special_function_cases(args.output)
    else: run(args.raw_dir, args.output, args.study_protocol)
