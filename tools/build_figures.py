#!/usr/bin/env python3
"""Build all current figures in a new isolated output directory."""
import argparse
import json
import sys
sys.dont_write_bytecode=True
from _figure_support import build_current

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,help='New directory under runs/ or outside the repository')
    args=parser.parse_args()
    try:
        result=build_current(args.output)
        print(json.dumps({'status':result['status'],'active_figures':result['active_figures'],
                          'companion_figures':result['companion_figures'],
                          'output':args.output,'repository_preserved':result['repository_integrity_after']['ok']},indent=2))
    except Exception as exc:
        print('Figure build failed: '+str(exc),file=sys.stderr)
        raise SystemExit(1)
