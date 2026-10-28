#!/usr/bin/env python3
"""Check production-summary invariance to changed diagnostic rule counts.

Only deep copies are perturbed. This intentionally calls the unchanged
production summarize routine; the independent statistics audit is separate.
The declared law, costs, and original input files are never modified.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from representation import load_representation
sys.path.insert(0, str(ROOT/'src'))
import joint_validation as production


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def object_sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=ROOT/'results/main')
    parser.add_argument('--output', type=Path, default=ROOT/'qa/production_no_refit.json')
    args = parser.parse_args()
    args.results = args.results.resolve()
    protocol_path = ROOT/'Protocol_J1_Replication.json'
    raw_path = args.results/'raw/J1_n25600.csv.gz'
    summary_path = args.results/'summary.json'
    receipt_path = ROOT/'integrity_record.json'
    receipt, current = load_representation(ROOT)
    immutable_paths = [protocol_path, raw_path, summary_path, receipt_path,
                       *[ROOT/relative for relative in current['source_sha256']]]
    before_files = {str(p.relative_to(ROOT)): sha(p) for p in immutable_paths}
    for relative, digest in current['source_sha256'].items():
        assert sha(ROOT/relative) == digest, relative
    assert sha(protocol_path) == current['protocol_sha256']
    protocol = json.loads(protocol_path.read_text())
    with gzip.open(raw_path, 'rt', newline='') as stream:
        rows = list(csv.DictReader(stream))
    integer_fields = ['n', 'replicate', 'residual_size', 'double_rows', 'G', 'G_prime',
                      'cyclic_components', 'largest_complex', 'cost_lower', 'cost_upper',
                      'solver_work']
    for row in rows:
        for field in integer_fields:
            row[field] = int(row[field])
    assert len(rows) == 20000 and [r['replicate'] for r in rows] == list(range(20000))
    rule_fields = ['rule_counts_d1', 'rule_counts_d3']
    baseline_vectors = [{field: json.loads(row[field]) for field in rule_fields} for row in rows]
    baseline_vectors_hash = object_sha(baseline_vectors)
    original_memory_hash = object_sha([rows, protocol])
    baseline = production.summarize('J1', 25600, rows, protocol)
    assert object_sha([rows, protocol]) == original_memory_hash
    assert baseline == json.loads(summary_path.read_text())[0]
    altered_vectors = deepcopy(baseline_vectors)
    altered_rows = deepcopy(rows)
    altered_protocol = deepcopy(protocol)
    changed_vectors = 0
    positive_codes = [0, 4, 5, 6]
    for original, altered, row in zip(baseline_vectors, altered_vectors, altered_rows):
        for field in rule_fields:
            original_counts = original[field]
            counts = altered[field]
            values = [counts[j] for j in positive_codes]
            for code, value in zip(positive_codes, values[1:] + values[:1]):
                counts[code] = value
            assert sum(counts) == sum(original_counts) == 12800
            assert all(counts[code] == 0 for code in [1, 2, 3, 7])
            changed_vectors += counts != original_counts
            row[field] = json.dumps(counts, separators=(',', ':'))
    assert changed_vectors == 40000
    for original, altered in zip(rows, altered_rows):
        assert {k: v for k, v in original.items() if k not in rule_fields} == {
            k: v for k, v in altered.items() if k not in rule_fields}
    assert altered_protocol == protocol
    altered_input_hash = object_sha([altered_rows, altered_protocol])
    perturbed = production.summarize('J1', 25600, altered_rows, altered_protocol)
    assert object_sha([altered_rows, altered_protocol]) == altered_input_hash
    assert object_sha([rows, protocol]) == original_memory_hash
    assert object_sha(baseline_vectors) == baseline_vectors_hash
    assert baseline == perturbed
    after_files = {str(p.relative_to(ROOT)): sha(p) for p in immutable_paths}
    assert before_files == after_files
    reference_keys = ['lambda', 'rho', 'contraction_coefficient', 'mu', 'reference_m2', 'reference_p0']
    report = dict(status='PASS', draws=20000, changed_diagnostic_vectors=changed_vectors,
                  perturbation='Rotate counts across the four positive-probability corner codes independently in both degree classes of every draw, preserving each vector total and all impossible categories; only deep copies are changed.',
                  checks=dict(recorded_production_summary_reproduced=True,
                              entire_production_summary_unchanged=True,
                              cost_summaries_unchanged=True, fixed_reference_unchanged=True,
                              declared_protocol_unchanged=True, all_nondiagnostic_fields_unchanged=True,
                              original_memory_inputs_unchanged=True, original_rule_vectors_unchanged=True,
                              perturbed_inputs_not_mutated_by_summarizer=True,
                              input_files_unchanged=True),
                  fixed_reference={key: baseline[key] for key in reference_keys},
                  original_rule_vectors_sha256=baseline_vectors_hash,
                  altered_rule_vectors_sha256=object_sha(altered_vectors),
                  production_summary_sha256=object_sha(baseline),
                  input_sha256=before_files,
                  test_source_sha256=sha(Path(__file__)),
                  scope='Dependency-separation regression check on the unchanged production summarizer; the separate independent statistics audit establishes the formulas and values.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': report['status'], 'output': str(args.output),
                      'changed_diagnostic_vectors': changed_vectors,
                      'input_files_unchanged': True}), flush=True)


if __name__ == '__main__':
    main()
