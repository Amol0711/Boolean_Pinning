"""Focused exact and seeded validation for general joint two-slot rule laws."""
from __future__ import annotations
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from fractions import Fraction
import gzip
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.stats import beta, norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'vendor'))
import engine
from scc_control import from_two_input_arrays, solve_network

FIELDS = ['family','n','replicate','network_sha256','residual_size','double_rows',
          'G','G_prime','certificate_cycle_count','cyclic_components','largest_complex',
          'cost_lower','cost_upper','exact_reason','stop_summary','solver_work','pins']


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def rational(x): return [x.numerator, x.denominator]
def corners(code): return [(code >> j) & 1 for j in range(3)]


def atomic_text(path, contents):
    temporary = path.with_name(path.name+'.crdownload')
    temporary.write_text(contents, encoding='utf-8')
    temporary.replace(path)


def write_json(path, value):
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True)+'\n')


def type_blocks(n, family):
    blocks = []
    for spec in family['types']:
        count = n*Fraction(*spec['fraction'])
        if count.denominator != 1: raise ValueError('Nonintegral type count')
        blocks.append((int(count),spec))
    if sum(c for c,_ in blocks)!=n or sum(c*t['degree'] for c,t in blocks)!=2*n:
        raise ValueError('Type counts or output stub balance invalid')
    return blocks


def generate(n, family, replicate, protocol):
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(
        [protocol['master_seed'],family['id'],n,replicate])))
    blocks = type_blocks(n,family)
    degrees = np.concatenate([np.full(c,t['degree'],dtype=np.int64) for c,t in blocks])
    parents = rng.permutation(np.repeat(np.arange(n,dtype=np.int64),degrees)).reshape(n,2)
    codes = np.concatenate([rng.choice(8,size=c,p=np.asarray(t['weights'],dtype=np.float64)/t['probability_denominator'])
                            for c,t in blocks])
    rules = ((codes[:,None] >> np.arange(3)) & 1).astype(np.uint8)
    return parents,rules,degrees


def analyze(family,n,replicate,parents,rules,degrees,caps):
    certified=engine.certify(parents,rules,degrees)
    R=np.flatnonzero(certified['R']).tolist(); components=None; mapping=None
    if certified['certified']:
        lower=upper=int(certified['cost']); pins=[int(v) for v in certified['pins']]
        Gp=True; cycles=lower; cyclic=lower; largest=0; work=0; stops={};reason='original_G'
    else:
        r=len(R); loc={v:i for i,v in enumerate(R)}
        p=[[loc.get(int(u),r) for u in parents[v]] for v in R]+[[r,r]]
        b=[rules[v].tolist() for v in R]+[[0,0,0]]
        solved=solve_network(from_two_input_arrays(p,b),support_mode='slot',**caps)
        assert solved['residual']==list(range(r)) and not solved['G']
        lower,upper=solved['lower'],solved['upper']
        pins=[R[v] for v in solved['feasible_pins']]
        Gp=solved['G_prime'];cycles=solved['certificate_cycle_count']
        cyclic=solved['cyclic_component_count'];work=solved['total_work']
        components=solved['components'];mapping=R+[None]
        largest=max([len(c['nodes']) for c in components if c['cyclic'] and not c['marked_slot_cycle']] or [0])
        stops=dict(Counter(c['reason'] for c in components if c['exact_cost'] is None))
        reason='weaker_G_prime' if Gp else ('SCC_exact' if lower==upper else 'bounded_unresolved')
    h=hashlib.sha256(np.asarray(parents,dtype='<i8').tobytes()+np.asarray(rules,dtype=np.uint8).tobytes()).hexdigest()
    pins=sorted(pins);assert len(pins)==upper
    row=dict(family=family,n=n,replicate=replicate,network_sha256=h,residual_size=len(R),
             double_rows=certified['double_rows'],G=int(certified['certified']),G_prime=int(Gp),
             certificate_cycle_count=cycles if Gp else '',cyclic_components=cyclic,largest_complex=largest,
             cost_lower=lower,cost_upper=upper,exact_reason=reason,
             stop_summary=json.dumps(stops,sort_keys=True,separators=(',',':')),solver_work=work,pins=';'.join(map(str,pins)))
    detail=dict(family=family,n=n,replicate=replicate,network_sha256=h,
                residual_rows=[[v,int(parents[v,0]),int(parents[v,1]),*[int(x) for x in rules[v]]] for v in R],
                local_to_original=mapping,components=components,G=bool(row['G']),G_prime=bool(Gp),
                certificate_cycle_count=cycles,lower=lower,upper=upper,pins=pins,exact_reason=reason)
    return row,detail


def stream(path):
    temporary=path.with_name(path.name+'.crdownload');raw=temporary.open('wb')
    gz=gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=6)
    return raw,io.TextIOWrapper(gz,encoding='utf-8',newline=''),temporary


def verify(path,count,kind):
    with gzip.open(path,'rt',encoding='utf-8',newline='') as f:
        rows=csv.DictReader(f) if kind=='raw' else (json.loads(line) for line in f)
        observed=0
        for observed,row in enumerate(rows,1):
            if int(row['replicate'])!=observed-1:raise RuntimeError('Record sequence mismatch')
        if observed!=count:raise RuntimeError('Incomplete stream')
    return digest(path)


def cp(k,N):
    return [0.0 if k==0 else float(beta.ppf(.025,k,N-k+1)),
            1.0 if k==N else float(beta.ppf(.975,k+1,N-k))]


def endpoint_moment(lo,hi,power):
    N=len(lo);x=[a**power for a in lo];y=[a**power for a in hi]
    ml=sum(x)/N;mu=sum(y)/N
    sl=math.sqrt(max(0,(sum(a*a for a in x)-N*ml*ml)/(N-1)/N))
    su=math.sqrt(max(0,(sum(a*a for a in y)-N*mu*mu)/(N-1)/N))
    z=float(norm.ppf(.975))
    return dict(lower=ml,upper=mu,se_lower=sl,se_upper=su,
                normal95_envelope=[max(0,ml-z*sl),mu+z*su])


def summarize(name,n,rows,protocol):
    N=len(rows);f=protocol['families'][name];lam=float(Fraction(*f['parameters_rational']['lambda']))
    rho=float(Fraction(*f['parameters_rational']['rho']));mu=-math.log1p(-lam)
    lo=[r['cost_lower'] for r in rows];hi=[r['cost_upper'] for r in rows]
    gf=sum(1-r['G'] for r in rows);gpf=sum(1-r['G_prime'] for r in rows)
    definite=sum(x==0 for x in hi);possible=sum(x==0 for x in lo)
    return dict(family=name,n=n,draws=N,**{'lambda':lam},rho=rho,contraction_coefficient=max(lam,rho),
                mu=mu,reference_m2=mu+mu*mu,reference_p0=1-lam,unresolved=sum(a!=b for a,b in zip(lo,hi)),
                G_failures=gf,G_prime_failures=gpf,G_failures_cp95=cp(gf,N),G_prime_failures_cp95=cp(gpf,N),
                mean=endpoint_moment(lo,hi,1),m2=endpoint_moment(lo,hi,2),
                P0=dict(definite_zero_count=definite,possible_zero_count=possible,lower=definite/N,upper=possible/N,
                        cp95_envelope=[cp(definite,N)[0],cp(possible,N)[1]]),
                max_residual=max(r['residual_size'] for r in rows),max_complex=max(r['largest_complex'] for r in rows),
                total_solver_work=sum(r['solver_work'] for r in rows),
                stop_reason_counts=dict(sum((Counter(json.loads(r['stop_summary'])) for r in rows),Counter())),
                exact_reason_counts=dict(Counter(r['exact_reason'] for r in rows)))


def cell(args):
    name,n,protocol,output=args;output=Path(output);N=protocol['replicates_per_cell']
    paths=[output/'raw'/f'{name}_n{n}.csv.gz',output/'details'/f'{name}_n{n}.jsonl.gz']
    opened=[stream(p) for p in paths];rows=[];started=time.perf_counter()
    try:
        writer=csv.DictWriter(opened[0][1],fieldnames=FIELDS);writer.writeheader()
        for rep in range(N):
            arrays=generate(n,protocol['families'][name],rep,protocol)
            row,detail=analyze(name,n,rep,*arrays,protocol['component_caps']);rows.append(row)
            writer.writerow(row);opened[1][1].write(json.dumps(detail,sort_keys=True,separators=(',',':'))+'\n')
    finally:
        for raw,text,_ in opened:text.close();raw.close()
    hashes={}
    for path,(_,_,temporary),kind in zip(paths,opened,['raw','details']):
        h=verify(temporary,N,kind);temporary.replace(path)
        if verify(path,N,kind)!=h:raise RuntimeError('Published stream changed')
        hashes[kind]=h
    summary=summarize(name,n,rows,protocol);summary['stream_sha256']=hashes
    return summary,dict(family=name,n=n,elapsed_seconds=time.perf_counter()-started)


def flatten(summary):
    result={k:v for k,v in summary.items() if not isinstance(v,(dict,list))}
    for label in ['mean','m2']:
        for k in ['lower','upper','se_lower','se_upper']:result[f'{label}_{k}']=summary[label][k]
        result[f'{label}_normal95_lo'],result[f'{label}_normal95_hi']=summary[label]['normal95_envelope']
    for k in ['definite_zero_count','possible_zero_count','lower','upper']:result[f'P0_{k}']=summary['P0'][k]
    result['P0_cp95_lo'],result['P0_cp95_hi']=summary['P0']['cp95_envelope']
    for k in ['G_failures','G_prime_failures']:result[k+'_cp95_lo'],result[k+'_cp95_hi']=summary[k+'_cp95']
    return result


def run_main(output,protocol,workers,gate):
    check=json.loads(gate.read_text())
    if check.get('status')!='PASS':raise ValueError('Independent exact validation gate must pass')
    for d in ['raw','details']:(output/d).mkdir(parents=True)
    engine.backend_info();started=time.perf_counter();summaries=[];resources=[]
    jobs=[(name,n,protocol,str(output)) for name in sorted(protocol['families']) for n in protocol['sizes']]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for future in as_completed([pool.submit(cell,job) for job in jobs]):
            summary,resource=future.result();summaries.append(summary);resources.append(resource)
            print(json.dumps({k:summary[k] for k in ['family','n','draws','G_failures','G_prime_failures','unresolved']}),flush=True)
    summaries.sort(key=lambda r:(r['family'],r['n']));resources.sort(key=lambda r:(r['family'],r['n']))
    integrity=[]
    for s in summaries:
        for kind,suffix in [('raw','csv.gz'),('details','jsonl.gz')]:
            p=output/kind/f"{s['family']}_n{s['n']}.{suffix}";h=verify(p,s['draws'],kind)
            if h!=s['stream_sha256'][kind]:raise RuntimeError('Completed stream changed after publication')
            integrity.append(dict(path=str(p.relative_to(output)),records=s['draws'],sha256=h))
    write_json(output/'summary.json',summaries)
    buf=io.StringIO(newline='');flat=[flatten(s) for s in summaries]
    writer=csv.DictWriter(buf,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    atomic_text(output/'summary.csv',buf.getvalue())
    write_json(output/'io_integrity.json',dict(status='PASS',streams=integrity))
    write_json(output/'resources.json',dict(workers=workers,total_elapsed_seconds=time.perf_counter()-started,cells=resources))
    write_json(output/'manifest.json',dict(status='completed',total_draws=sum(s['draws'] for s in summaries),
               protocol_sha256=digest(ROOT/'Protocol_JointRules.json'),independent_gate_sha256=digest(gate),
               source_sha256={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),ROOT/'vendor/engine.py',ROOT/'vendor/engine_core.cpp',ROOT/'vendor/scc_control.py']},
               backend=engine.backend_info(),numpy=np.__version__,unresolved=sum(s['unresolved'] for s in summaries)))


def run_exact(output,protocol):
    fixtures=dict(protocol['families'])
    for name,weights in [('ZERO_AND',[40,0,0,0,60,0,0,0]),('ZERO',[100,0,0,0,0,0,0,0]),('ALL_AND',[0,0,0,0,100,0,0,0])]:
        fixtures[name]={'types':[{'degree':2,'fraction':[1,1],'weights':weights,'probability_denominator':100}]}
    outcomes=[];summaries=[];support=[]
    for code in range(8):
        table=[0]+corners(code)
        values=[int(any(table[e] for e in range(4) if e & ~mask==0)) for mask in range(4)]
        support.append(dict(corner_code=code,corners=corners(code),active_mask_support=values))
    for name,family in fixtures.items():
        specs=[t for c,t in type_blocks(2,family) for _ in range(c)]
        ds=[t['degree'] for t in specs];sources=[i for i,d in enumerate(ds) for _ in range(d)]
        matchings=sorted(set(itertools.permutations(sources)))
        choices=[[(c,Fraction(w,t['probability_denominator'])) for c,w in enumerate(t['weights']) if w] for t in specs]
        sums={k:Fraction(0) for k in ['weight','C1','C2','C1_factorial2','cost','G','G_prime']};dist=Counter();records=0
        for matching in matchings:
            parents=[list(matching[:2]),list(matching[2:])]
            for selected in itertools.product(*choices):
                codes=[v[0] for v in selected];weight=selected[0][1]*selected[1][1]/len(matchings)
                rules=[corners(c) for c in codes]
                solved=solve_network(from_two_input_arrays(parents,rules),support_mode='slot',**protocol['component_caps'])
                assert solved['exact_cost'] is not None
                c1=sum(rules[v][slot] for v in range(2) for slot in range(2) if parents[v][slot]==v)
                c2=sum(rules[0][s]*rules[1][t] for s in range(2) for t in range(2) if parents[0][s]==1 and parents[1][t]==0)
                record=dict(family=name,n=2,parents=parents,corner_codes=codes,weight_rational=rational(weight),
                            C1=c1,C2=c2,cost=solved['exact_cost'],feasible_pins=solved['feasible_pins'],G=solved['G'],G_prime=solved['G_prime'],residual=solved['residual'])
                if solved['G_prime']:assert solved['exact_cost']==c1+c2
                outcomes.append(record);records+=1;dist[solved['exact_cost']]+=weight
                for key,value in {'weight':1,'C1':c1,'C2':c2,'C1_factorial2':c1*(c1-1),'cost':solved['exact_cost'],'G':int(solved['G']),'G_prime':int(solved['G_prime'])}.items():sums[key]+=weight*value
        a=[sum(p*((c&1)+((c>>1)&1)) for c,p in ch) for ch in choices]
        q12=[sum(p for c,p in ch if c&3==3) for ch in choices]
        cross=sum(ds[v]*a[v]*ds[w]*a[w] for v in range(2) for w in range(2) if v!=w)
        expected={'C1':sum(ds[v]*a[v] for v in range(2))/4,'C2':cross/24,
                  'C1_factorial2':(cross+2*sum(ds[v]*(ds[v]-1)*q12[v] for v in range(2)))/12}
        assert sums['weight']==1 and all(sums[k]==v for k,v in expected.items())
        summaries.append(dict(family=name,realizations=records,expected_rational={k:rational(v) for k,v in expected.items()},
                              enumerated_rational={k:rational(v) for k,v in sums.items()},cost_distribution_rational={str(k):rational(v) for k,v in sorted(dist.items())}))
    atomic_text(output/'outcomes.jsonl',''.join(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n' for r in outcomes))
    write_json(output/'support_patterns.json',support)
    write_json(output/'summary.json',dict(status='PASS',protocol_sha256=digest(ROOT/'Protocol_JointRules.json'),realizations=len(outcomes),fixtures=summaries,source_sha256=digest(Path(__file__))))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['exact','main']);parser.add_argument('--output',type=Path)
    parser.add_argument('--workers',type=int,default=4);parser.add_argument('--gate',type=Path)
    args=parser.parse_args()
    if args.workers<1:parser.error('workers must be positive')
    if args.mode=='main' and args.gate is None:parser.error('main requires --gate independent-PASS.json')
    output=args.output or ROOT/'results'/args.mode
    if output.exists():parser.error('Output exists; refusing overwrite')
    protocol=json.loads((ROOT/'Protocol_JointRules.json').read_text())
    if protocol['fixed_design'] is not True:raise ValueError('Protocol must use a fixed design')
    output.mkdir(parents=True)
    if args.mode=='exact':run_exact(output,protocol)
    else:run_main(output,protocol,args.workers,args.gate)


if __name__=='__main__':main()
