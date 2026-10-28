#!/usr/bin/env python3
"""Check the current computational-figure inputs in a temporary copy."""
import json
import sys
sys.dont_write_bytecode=True
from _figure_support import verify_current

if __name__=='__main__':
    try:
        if len(sys.argv)!=1:
            raise ValueError('Usage: python tools/verify_figures.py')
        result=verify_current()
        print(json.dumps(result,indent=2,sort_keys=True))
    except Exception as exc:
        print('Figure verification failed: '+str(exc),file=sys.stderr)
        raise SystemExit(1)
