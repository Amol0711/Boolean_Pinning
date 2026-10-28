"""J1 replication with deterministic random streams and complete sample accounting."""
from __future__ import annotations
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import io
import itertools
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import scipy
import joint_validation as jv
from scc_control import from_two_input_arrays, _sccs, _cycle_info, solve_network

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from representation import load_representation, protocol_reference
PROTOCOL=ROOT/'Protocol_J1_Replication.json'
FIELDS=jv.FIELDS+['seed_entropy','rule_counts_d1','rule_counts_d3','component_sizes','cyclic_component_sizes']
SOURCES=['src/run_replication.py','src/joint_validation.py','vendor/engine.py','vendor/engine_core.cpp','vendor/scc_control.py']

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read_json(path): return json.loads(Path(path).read_text())
def compact(value): return json.dumps(value,separators=(',',':'),sort_keys=True)

def deterministic_preflight():
    sys.path.insert(0,str(ROOT/'verification'))
    from independent_network_oracle import exact_cost,witness_success
    checked=0
    for parents in itertools.product(range(2),repeat=4):
        p=np.array(parents,dtype=np.int64).reshape(2,2)
        degrees=np.bincount(p.ravel(),minlength=2)
        for codes in itertools.product(range(8),repeat=2):
            b=((np.array(codes)[:,None] >> np.arange(3)) & 1).astype(np.uint8)
            fast=jv.engine.certify(p,b,degrees,backend='cpp')
            slow=jv.engine.certify(p,b,degrees,backend='python')
            for key in ['R','pins','residual_parent']: assert np.array_equal(fast[key],slow[key]),(parents,codes,key)
            for key in ['double_rows','certified','cost']: assert fast[key]==slow[key]
            rows=[[v,int(p[v,0]),int(p[v,1]),*map(int,b[v])] for v in range(2)]
            oracle=exact_cost(rows,[0,1])
            produced=solve_network(from_two_input_arrays(p,b),max_component_size=18,max_pins=18,max_work=10000000)
            assert produced['exact_cost']==oracle['cost']
            assert witness_success(rows,[0,1],produced['feasible_pins'])['success']
            checked+=1
    original=read_json(ROOT/'original/Protocol_JointRules.json')
    with gzip.open(ROOT/'original/raw/J1_n25600.csv.gz','rt') as f:
        old=list(csv.DictReader(f))
    for rep in [0,999,1999]:
        p,b,d=jv.generate(25600,original['families']['J1'],rep,original)
        row,_=jv.analyze('J1',25600,rep,p,b,d,original['component_caps'])
        assert row['network_sha256']==old[rep]['network_sha256']
        for key in ['cost_lower','cost_upper','G','G_prime']: assert row[key]==int(old[rep][key])
    return {'status':'PASS','deterministic_two_node_networks':checked,'original_draws_replayed':[0,999,1999],
            'production_seed_draws':0,'backend':jv.engine.backend_info()}

def verify_integrity():
    receipt,current=load_representation(ROOT)
    assert np.__version__==receipt['software']['numpy']
    assert scipy.__version__==receipt['software']['scipy']
    return receipt,current

def one_draw(rep,protocol):
    n=25600
    p,b,d=jv.generate(n,protocol['families']['J1'],rep,protocol)
    row,detail=jv.analyze('J1',n,rep,p,b,d,protocol['component_caps'])
    codes=b[:,0]+2*b[:,1]+4*b[:,2]
    counts=[np.bincount(codes[:12800],minlength=8).tolist(),np.bincount(codes[12800:],minlength=8).tolist()]
    R=[x[0] for x in detail['residual_rows']];r=len(R);loc={v:i for i,v in enumerate(R)}
    localp=[[loc.get(int(u),r) for u in p[v]] for v in R]+[[r,r]]
    localb=[b[v].tolist() for v in R]+[[0,0,0]]
    net=from_two_input_arrays(localp,localb)
    components=_sccs(net,list(range(r)))
    sizes=sorted(map(len,components))
    cyclic_sizes=sorted(len(c) for c in components if _cycle_info(net,c)['cyclic'])
    extras={'seed_entropy':[protocol['master_seed'],protocol['families']['J1']['id'],n,rep],
            'rule_counts_d1':counts[0],'rule_counts_d3':counts[1],
            'component_sizes':sizes,'cyclic_component_sizes':cyclic_sizes}
    assert sum(sizes)==r
    for count in counts: assert sum(count)==12800 and all(count[c]==0 for c in [1,2,3,7])
    detail.update(extras);row.update({key:compact(value) for key,value in extras.items()})
    return row,detail

def shard(args):
    first,last,output=args;output=Path(output);protocol=read_json(PROTOCOL)
    verify_integrity();start=time.perf_counter()
    paths=[output/f'{first:05d}_{last:05d}.csv.gz',output/f'{first:05d}_{last:05d}.jsonl.gz']
    handles=[jv.stream(p) for p in paths]
    try:
        writer=csv.DictWriter(handles[0][1],fieldnames=FIELDS);writer.writeheader()
        for rep in range(first,last):
            row,detail=one_draw(rep,protocol)
            writer.writerow(row);handles[1][1].write(compact(detail)+'\n')
    finally:
        for raw,text,_ in handles:text.close();raw.close()
    for path,(_,_,temporary),kind in zip(paths,handles,['raw','detail']):
        with gzip.open(temporary,'rt') as f:
            rows=csv.DictReader(f) if kind=='raw' else (json.loads(line) for line in f)
            assert [int(row['replicate']) for row in rows]==list(range(first,last))
        temporary.replace(path)
    return {'first':first,'last_exclusive':last,'draws':last-first,'elapsed_seconds':time.perf_counter()-start,
            'raw_sha256':sha(paths[0]),'detail_sha256':sha(paths[1])}

def run(workers,output):
    receipt,current=verify_integrity();protocol=read_json(PROTOCOL)
    if output.exists(): raise RuntimeError('Output exists; refusing overwrite or outcome-based continuation')
    for directory in ['raw','details','shards']:(output/directory).mkdir(parents=True,exist_ok=True)
    jv.write_json(output/'run_started.json',{'started_at_utc':datetime.now(timezone.utc).isoformat(),
        'integrity_record_sha256':sha(ROOT/'integrity_record.json'),'workers':workers,'requested_draws':20000})
    started=time.perf_counter();completed=[]
    jobs=[(first,min(20000,first+1000),str(output/'shards')) for first in range(0,20000,1000)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for future in as_completed([pool.submit(shard,job) for job in jobs]):
            completed.append(future.result())
            print(json.dumps({'completed_draws':sum(x['draws'] for x in completed),'requested_draws':20000}),flush=True)
    completed.sort(key=lambda x:x['first'])
    paths=[output/'raw/J1_n25600.csv.gz',output/'details/J1_n25600.jsonl.gz']
    handles=[jv.stream(p) for p in paths];raw_rows=[]
    try:
        writer=csv.DictWriter(handles[0][1],fieldnames=FIELDS);writer.writeheader()
        for item in completed:
            stem=f"{item['first']:05d}_{item['last_exclusive']:05d}"
            with gzip.open(output/'shards'/f'{stem}.csv.gz','rt') as stream:
                for row in csv.DictReader(stream):
                    writer.writerow(row)
                    for key in ['replicate','n','residual_size','double_rows','G','G_prime','cyclic_components','largest_complex','cost_lower','cost_upper','solver_work']:row[key]=int(row[key])
                    raw_rows.append(row)
            with gzip.open(output/'shards'/f'{stem}.jsonl.gz','rt') as stream:
                for line in stream:handles[1][1].write(line)
    finally:
        for raw,text,_ in handles:text.close();raw.close()
    integrity=[]
    for path,(_,_,temporary),kind in zip(paths,handles,['raw','details']):
        jv.verify(temporary,20000,kind);temporary.replace(path)
        integrity.append({'path':str(path.relative_to(output)),'records':20000,'sha256':jv.verify(path,20000,kind)})
    assert len(raw_rows)==20000 and [r['replicate'] for r in raw_rows]==list(range(20000))
    summary=jv.summarize('J1',25600,raw_rows,protocol)
    jv.write_json(output/'summary.json',[summary])
    with (output/'summary.csv').open('w',newline='') as f:
        flat=jv.flatten(summary);writer=csv.DictWriter(f,fieldnames=list(flat));writer.writeheader();writer.writerow(flat)
    counts={str(d):[sum(json.loads(r[f'rule_counts_d{d}'])[c] for r in raw_rows) for c in range(8)] for d in [1,3]}
    jv.write_json(output/'rule_counts.json',{'order':'corner_code0..7','draws':20000,'rows_per_degree':256000000,'counts':counts})
    jv.write_json(output/'io_integrity.json',{'status':'PASS','streams':integrity})
    jv.write_json(output/'resources.json',{'workers':workers,'elapsed_seconds':time.perf_counter()-started,'shards':completed})
    verify_integrity()
    jv.write_json(output/'manifest.json',{'status':'completed','total_draws':20000,'unresolved':summary['unresolved'],
        'protocol_sha256':sha(PROTOCOL),'integrity_record_sha256':sha(ROOT/'integrity_record.json'),
        'source_sha256':current['source_sha256'],'software':receipt['software'],
        'independent_gate_sha256':current['independent_preflight_sha256'],'backend':jv.engine.backend_info(),
        'protocol_record':protocol_reference(receipt),
        'integrity_note':'Current executable and protocol hashes are checked before and after the run.'})
    print(json.dumps({'status':'COMPLETED','draws':20000,'unresolved':summary['unresolved'],'elapsed_seconds':time.perf_counter()-started}),flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['run'])
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--output',type=Path,default=ROOT/'results/main')
    args=parser.parse_args()
    if not 1<=args.workers<=8:parser.error('workers must be 1..8')
    run(args.workers,args.output)

if __name__=='__main__':main()
