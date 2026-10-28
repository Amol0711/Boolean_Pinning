"""Repeat the recorded exact/main JointRules computations and compare scientific payloads."""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--exact-only',action='store_true')
    parser.add_argument('--gate',type=Path,default=ROOT.parent/'qa/independent_exact_checks.json')
    args=parser.parse_args()
    if args.workers<1:parser.error('workers must be positive')
    if args.out is not None and args.out.exists():parser.error('--out must not exist')
    out=args.out.resolve() if args.out else Path(tempfile.mkdtemp(prefix='joint_rules_reproduction_'))/'replay'
    out.mkdir(parents=True)
    modes=['exact'] if args.exact_only else ['exact','main']
    comparisons=[]
    for mode in modes:
        command=[sys.executable,str(ROOT/'src/joint_validation.py'),mode,'--output',str(out/mode),'--workers',str(args.workers)]
        if mode=='main':command+=['--gate',str(args.gate.resolve())]
        subprocess.run(command,check=True)
        relative=([Path(name) for name in ['outcomes.jsonl','summary.json','support_patterns.json']]
                  if mode=='exact' else [Path(name) for name in ['summary.json','summary.csv']]+
                  sorted(p.relative_to(ROOT/'results/main') for kind in ['raw','details'] for p in (ROOT/'results/main'/kind).glob('*.gz')))
        if mode=='main' and len(relative)!=20:raise RuntimeError('Expected18 streams plus2 summaries')
        for path in relative:
            expected=(ROOT/'results'/mode/path).read_bytes();actual=(out/mode/path).read_bytes()
            compressed=str(path).endswith('.gz')
            ep=gzip.decompress(expected) if compressed else expected
            ap=gzip.decompress(actual) if compressed else actual
            comparisons.append(dict(path=str(Path(mode)/path),byte_identical=expected==actual,
                                    payload_identical=ep==ap,reference_sha256=hashlib.sha256(expected).hexdigest(),replay_sha256=hashlib.sha256(actual).hexdigest()))
    passed=all(r['byte_identical'] and r['payload_identical'] for r in comparisons)
    report=dict(status='PASS' if passed else 'FAIL',exact_only=args.exact_only,scientific_files=len(comparisons),files=comparisons,
                timing_and_manifest_metadata_compared=False)
    temporary=out/'reproduction_comparison.json.crdownload';temporary.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    temporary.replace(out/'reproduction_comparison.json')
    print(json.dumps(dict(status=report['status'],scientific_files=len(comparisons),report=str(out/'reproduction_comparison.json'))),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
