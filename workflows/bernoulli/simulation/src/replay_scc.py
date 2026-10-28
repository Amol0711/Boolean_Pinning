"""Original-seed configuration-network replay with SCC-localized bounds.

All scientific output ordering is deterministic; run resources are separate.
Baseline files and their generated scientific records are never overwritten.
"""
from __future__ import annotations
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'vendor'))
import engine
from scc_control import from_two_input_arrays, solve_network

FIELDS = ['mixture', 'n', 'replicate', 'network_sha256', 'residual_size', 'double_rows',
          'G', 'G_prime', 'certificate_cycle_count', 'cyclic_components', 'largest_complex',
          'scc_lower', 'scc_upper', 'original_lower', 'original_upper', 'cost_lower', 'cost_upper',
          'exact_reason', 'stop_summary', 'solver_work', 'pins', 'bound_witness_source']


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def network_digest(parents, rules):
    h = hashlib.sha256()
    h.update(np.asarray(parents, dtype='<i8').tobytes(order='C'))
    h.update(np.asarray(rules, dtype=np.uint8).tobytes(order='C'))
    return h.hexdigest()


def original_records(mixture, n):
    path = ROOT / 'baseline/results/refined/raw' / f'{mixture}_n{n}.csv.gz'
    with gzip.open(path, 'rt', newline='') as f:
        return list(csv.DictReader(f))


def analyze_draw(mixture, n, original, caps):
    rep = int(original['replicate'])
    rng = engine.make_rng(2026101301, engine.MIXTURE_IDS[mixture], n, rep)
    parents, rules, degrees, _ = engine.generate_network(n, mixture, rng)
    certified = engine.certify(parents, rules, degrees)
    for key in ('residual_size', 'double_rows'):
        assert certified[key] == int(original[key]), (mixture, n, rep, key)
    assert int(certified['certified']) == int(original['certified'])
    R = np.flatnonzero(certified['R']).tolist()
    components = None
    local_to_original = R + [None]
    oldlo, oldhi = int(original['cost_lower']), int(original['cost_upper'])
    oldpins = [] if not original['pins'] else [int(x) for x in original['pins'].split(';')]
    assert len(oldpins) == oldhi
    if certified['certified']:
        slo = shi = int(certified['cost'])
        sp = [int(x) for x in certified['pins']]
        assert slo == oldlo == oldhi
        Gprime = True; cyclecount = slo; cyclic = slo; largest = 0; work = 0
        stops = {}
    else:
        # Exterior nominal slots point to one new constant-zero coordinate.
        # Both slots are retained, including equal physical source labels.
        r = len(R)
        loc = {v: j for j, v in enumerate(R)}
        p = [[loc.get(int(u), r) for u in parents[v]] for v in R] + [[r, r]]
        b = [rules[v].tolist() for v in R] + [[0, 0, 0]]
        net = from_two_input_arrays(p, b)
        solved = solve_network(net, support_mode='slot', **caps)
        assert solved['residual'] == list(range(r))
        assert not solved['G']
        slo, shi = solved['lower'], solved['upper']
        sp = [R[i] for i in solved['feasible_pins']]
        Gprime = solved['G_prime']; cyclecount = solved['certificate_cycle_count']
        cyclic = solved['cyclic_component_count']; work = solved['total_work']
        components = solved['components']
        largest = max([len(c['nodes']) for c in components if c['cyclic'] and not c['marked_slot_cycle']] or [0])
        stops = dict(Counter(c['reason'] for c in components if c['exact_cost'] is None))
    lower, upper = max(oldlo, slo), min(oldhi, shi)
    assert lower <= upper, (mixture, n, rep, oldlo, oldhi, slo, shi)
    if oldlo == oldhi: assert lower == upper == oldlo
    if oldhi < shi:
        pins, source = oldpins, 'baseline_feasible_witness'
    else:
        pins, source = sorted(sp), 'SCC_or_G_feasible_witness'
    assert len(pins) == upper
    if lower == upper:
        reason = 'original_G' if certified['certified'] else ('weaker_G_prime' if Gprime else
                 ('SCC_exact' if slo == shi else ('preserved_original_exact' if oldlo == oldhi else 'intersected_bounds_exact')))
    else: reason = 'bounded_unresolved'
    row = dict(mixture=mixture, n=n, replicate=rep, network_sha256=network_digest(parents, rules),
               residual_size=len(R), double_rows=certified['double_rows'], G=int(certified['certified']),
               G_prime=int(Gprime), certificate_cycle_count=cyclecount if Gprime else '',
               cyclic_components=cyclic, largest_complex=largest, scc_lower=slo, scc_upper=shi,
               original_lower=oldlo, original_upper=oldhi, cost_lower=lower, cost_upper=upper,
               exact_reason=reason, stop_summary=json.dumps(stops, sort_keys=True, separators=(',', ':')),
               solver_work=work, pins=';'.join(map(str, pins)), bound_witness_source=source)
    detail = {'mixture': mixture, 'n': n, 'replicate': rep, 'network_sha256': row['network_sha256'],
              'residual_rows': [[v, int(parents[v, 0]), int(parents[v, 1]), *[int(x) for x in rules[v]]] for v in R],
              'local_to_original': local_to_original if components is not None else None,
              'components': components, 'original_record': original,
              'G': bool(certified['certified']), 'G_prime': Gprime,
              'certificate_cycle_count': cyclecount, 'lower': lower, 'upper': upper,
              'pins': pins, 'bound_witness_source': source, 'exact_reason': reason}
    return row, detail


def write_stream(path):
    # The shared-workspace synchronizer excludes *.crdownload. Keeping an
    # open stream out of synchronization also protects its pathname from
    # replacement by a partial remote snapshot while its descriptor is open.
    temporary = path.with_name(path.name + '.crdownload')
    raw = temporary.open('wb')
    gz = gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=6)
    return raw, io.TextIOWrapper(gz, encoding='utf-8', newline=''), temporary


def verify_stream(path, count, kind):
    """Read through the gzip trailer/CRC and check every record identity."""
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as f:
        rows = csv.DictReader(f) if kind == 'raw' else (json.loads(line) for line in f)
        seen = 0
        for seen, row in enumerate(rows, 1):
            if int(row['replicate']) != seen - 1:
                raise RuntimeError(f'Nonconsecutive record in {path}')
        if seen != count:
            raise RuntimeError(f'Incomplete {path}: {seen} records, expected {count}')
    return digest(path)


def atomic_text(path, contents):
    temporary = path.with_name(path.name + '.crdownload')
    with temporary.open('w', encoding='utf-8', newline='') as f:
        f.write(contents)
        f.flush()
        os.fsync(f.fileno())
    temporary.replace(path)


def process_cell(args):
    mixture, n, count, caps, output = args
    output = Path(output)
    original = original_records(mixture, n)
    assert len(original) == 10000
    assert all(row['mixture'] == mixture and int(row['n']) == n
               and int(row['replicate']) == i for i, row in enumerate(original))
    selected = original[:count]
    summary = {'mixture': mixture, 'n': n, 'draws': count, 'G_failures': 0,
               'G_prime_failures': 0, 'unresolved': 0, 'original_unresolved': 0,
               'newly_resolved': 0, 'max_residual': 0, 'max_complex': 0,
               'total_solver_work': 0, 'stop_reason_counts': {}, 'exact_reason_counts': {}}
    stops, reasons = Counter(), Counter()
    raw_path = output / 'raw' / f'{mixture}_n{n}.csv.gz'
    detail_path = output / 'details' / f'{mixture}_n{n}.jsonl.gz'
    rawfile, rawtext, raw_temporary = write_stream(raw_path)
    detailfile, detailtext, detail_temporary = write_stream(detail_path)
    started = time.perf_counter()
    try:
        writer = csv.DictWriter(rawtext, fieldnames=FIELDS); writer.writeheader()
        for old in selected:
            row, detail = analyze_draw(mixture, n, old, caps)
            writer.writerow(row)
            detailtext.write(json.dumps(detail, sort_keys=True, separators=(',', ':')) + '\n')
            summary['G_failures'] += 1-row['G']
            summary['G_prime_failures'] += 1-row['G_prime']
            summary['unresolved'] += row['cost_lower'] != row['cost_upper']
            old_unresolved = row['original_lower'] != row['original_upper']
            summary['original_unresolved'] += old_unresolved
            summary['newly_resolved'] += old_unresolved and row['cost_lower'] == row['cost_upper']
            summary['max_residual'] = max(summary['max_residual'], row['residual_size'])
            summary['max_complex'] = max(summary['max_complex'], row['largest_complex'])
            summary['total_solver_work'] += row['solver_work']
            stops.update(json.loads(row['stop_summary'])); reasons[row['exact_reason']] += 1
    finally:
        rawtext.close(); rawfile.close(); detailtext.close(); detailfile.close()
    # Never publish a partial gzip file. Validate before and after atomic
    # publication, then validate every stream again before run completion.
    raw_hash = verify_stream(raw_temporary, count, 'raw')
    detail_hash = verify_stream(detail_temporary, count, 'details')
    raw_temporary.replace(raw_path); detail_temporary.replace(detail_path)
    if verify_stream(raw_path, count, 'raw') != raw_hash or verify_stream(detail_path, count, 'details') != detail_hash:
        raise RuntimeError('Published stream differs from verified temporary file')
    elapsed = time.perf_counter()-started
    summary['stop_reason_counts'] = dict(sorted(stops.items()))
    summary['exact_reason_counts'] = dict(sorted(reasons.items()))
    summary['raw_sha256'] = raw_hash; summary['detail_sha256'] = detail_hash
    return summary, {'mixture': mixture, 'n': n, 'draws': count, 'elapsed_seconds': elapsed}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pilot', action='store_true')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--work-cap', type=int, default=10_000_000)
    args = parser.parse_args()
    if args.workers < 1 or args.work_cap < 0: parser.error('invalid worker count or work cap')
    caps = {'max_component_size': 18, 'max_pins': 18, 'max_work': args.work_cap}
    output = args.output or ROOT/'results'/('pilot' if args.pilot else 'main')
    if output.exists(): parser.error('Output directory exists; refusing overwrite')
    output.mkdir(parents=True)
    (output/'raw').mkdir(); (output/'details').mkdir()
    protocol = json.loads((ROOT/'baseline/Protocol_Baseline.json').read_text())
    count = 100 if args.pilot else 10000
    if not args.pilot:
        declared = json.loads((ROOT/'Protocol_Bernoulli.json').read_text())
        assert declared['component_caps'] == caps and declared['replicates_per_cell'] == count
        assert declared['seed_scheme'] == protocol['seed_scheme']
        assert declared['frozen_before_main_replay'] is True
    engine.backend_info()  # compile before multiprocessing, avoiding simultaneous writes
    resources, summaries = [], []
    start = time.perf_counter()
    jobs = [(m, n, count, caps, str(output)) for m in 'ABCDE' for n in protocol['sizes']]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        pending = {pool.submit(process_cell, job): job[:2] for job in jobs}
        for future in as_completed(pending):
            summary, resource = future.result()
            summaries.append(summary); resources.append(resource)
            print(json.dumps({'mixture':summary['mixture'],'n':summary['n'],'draws':count,
                              'G_failures':summary['G_failures'],'G_prime_failures':summary['G_prime_failures'],
                              'unresolved':summary['unresolved'],'seconds':round(resource['elapsed_seconds'],2)}),flush=True)
    summaries.sort(key=lambda r:(r['mixture'],r['n']));resources.sort(key=lambda r:(r['mixture'],r['n']))
    integrity = []
    for summary in summaries:
        name = f"{summary['mixture']}_n{summary['n']}"
        for kind, suffix, field in [('raw','csv.gz','raw_sha256'),('details','jsonl.gz','detail_sha256')]:
            path = output/kind/f'{name}.{suffix}'
            sha = verify_stream(path, count, kind)
            if sha != summary[field]:
                raise RuntimeError(f'Completed stream changed after publication: {path}')
            integrity.append({'path':str(path.relative_to(output)), 'records':count, 'sha256':sha})
    atomic_text(output/'io_integrity.json', json.dumps({'status':'PASS','streams':integrity},indent=2,sort_keys=True)+'\n')
    atomic_text(output/'operational_summary.json', json.dumps(summaries, indent=2, sort_keys=True)+'\n')
    csv_buffer=io.StringIO(newline='')
    writer=csv.DictWriter(csv_buffer,fieldnames=list(summaries[0]));writer.writeheader()
    for row in summaries:
        writer.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,dict) else v for k,v in row.items()})
    atomic_text(output/'operational_summary.csv',csv_buffer.getvalue())
    meta={'status':'completed','pilot':args.pilot,'total_draws':sum(r['draws'] for r in summaries),
          'component_caps':caps,'engine':engine.backend_info(),
          'source_hashes':{str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),ROOT/'vendor/engine.py',ROOT/'vendor/engine_core.cpp',ROOT/'vendor/scc_control.py']},
          'original_protocol_sha256':digest(ROOT/'baseline/Protocol_Baseline.json'),
          'Bernoulli_protocol_sha256': None if args.pilot else digest(ROOT/'Protocol_Bernoulli.json'),
          'original_refined_raw_sha256': {
              str(p.relative_to(ROOT)): digest(p)
              for p in sorted((ROOT/'baseline/results/refined/raw').glob('*.csv.gz'))},
          'original_G_failure_count':sum(r['G_failures'] for r in summaries),
          'weaker_G_prime_failure_count':sum(r['G_prime_failures'] for r in summaries),
          'unresolved_count':sum(r['unresolved'] for r in summaries)}
    atomic_text(output/'resources.json',json.dumps({'total_elapsed_seconds':time.perf_counter()-start,'workers':args.workers,'cells':resources},indent=2,sort_keys=True)+'\n')
    atomic_text(output/'run_manifest.json',json.dumps(meta,indent=2,sort_keys=True)+'\n')
    print(json.dumps(meta),flush=True)


if __name__ == '__main__': main()
