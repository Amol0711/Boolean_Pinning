"""Check every production tiny joint-law outcome without production imports."""
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import itertools as it
import json
from pathlib import Path
import numpy as np
from independent_exact_checks import cycles, code_bits
from independent_network_oracle import (exact_cost, export_rows, graph_certificate,
    literal_support, physical_successors, succeeds)

ROOT=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    directory=ROOT/'validation/results/exact'
    protocol=json.loads((ROOT/'validation/Protocol_JointRules.json').read_text())
    fixtures=dict(protocol['families'])
    for name,law in [('ZERO_AND',{0:40,4:60}),('ZERO',{0:100}),('ALL_AND',{4:100})]:
        fixtures[name]={'types':[{'degree':2,'fraction':[1,1],
            'weights':[law.get(c,0) for c in range(8)],'probability_denominator':100}]}
    expected={}
    for name,family in fixtures.items():
        specs=[]
        for row in family['types']:
            count=2*Fraction(*row['fraction']); assert count.denominator==1
            specs.extend([row]*int(count))
        source=[v for v,row in enumerate(specs) for _ in range(row['degree'])]
        labeled=Counter(tuple(source[j] for j in perm) for perm in it.permutations(range(4)))
        for flat,multiplicity in labeled.items():
            for codes in it.product(*[[c for c,w in enumerate(row['weights']) if w] for row in specs]):
                weight=Fraction(multiplicity,24)
                for c,row in zip(codes,specs):
                    weight*=Fraction(row['weights'][c],row['probability_denominator'])
                expected[(name,flat,codes)]=weight
    seen=set();counts=Counter();sums=defaultdict(lambda:defaultdict(Fraction));dists=defaultdict(Counter)
    for line in (directory/'outcomes.jsonl').read_text().splitlines():
        row=json.loads(line);name=row['family'];n=row['n'];assert n==2
        key=(name,tuple(x for pair in row['parents'] for x in pair),tuple(row['corner_codes']))
        assert key not in seen;seen.add(key)
        weight=Fraction(*row['weight_rational']);assert weight==expected[key]
        parents=np.array(row['parents'],dtype=np.int64)
        corners=np.array([code_bits(c) for c in row['corner_codes']],dtype=np.uint8)
        residual,_=literal_support(parents,corners)
        all_rows=export_rows(parents,corners,list(range(n)))
        residual_rows=export_rows(parents,corners,residual)
        certificate=graph_certificate(residual_rows)
        assert row['residual']==residual
        assert row['G']==certificate['G'] and row['G_prime']==certificate['G_prime']
        oracle=exact_cost(all_rows,list(range(n)))
        assert oracle['complete'] and oracle['cost']==row['cost']
        pins=row['feasible_pins'];assert len(pins)==len(set(pins))==row['cost']
        assert succeeds(physical_successors(all_rows,list(range(n))),sum(1 << v for v in pins))
        marked=Counter(len(c[0]) for c in cycles(row['parents'],row['corner_codes']))
        assert row['C1']==marked[1] and row['C2']==marked[2]
        if name=='ALL_AND':
            assert row['cost']>0 and row['C1']==row['C2']==0
        for k,v in {'weight':1,'C1':marked[1],'C2':marked[2],
            'C1_factorial2':marked[1]*(marked[1]-1),'cost':row['cost'],
            'G':int(certificate['G']),'G_prime':int(certificate['G_prime'])}.items():
            sums[name][k]+=weight*v
        dists[name][str(row['cost'])]+=weight;counts[name]+=1
    assert seen==set(expected)
    summary=json.loads((directory/'summary.json').read_text())
    assert summary['status']=='PASS' and summary['realizations']==len(seen)
    for row in summary['fixtures']:
        name=row['family'];assert row['realizations']==counts[name]
        assert sums[name]['weight']==1
        for k,v in row['enumerated_rational'].items():assert sums[name][k]==Fraction(*v)
        for k,v in row['expected_rational'].items():assert sums[name][k]==Fraction(*v)
        assert dists[name]=={k:Fraction(*v) for k,v in row['cost_distribution_rational'].items()}
    patterns=json.loads((directory/'support_patterns.json').read_text())
    assert sorted(x['corner_code'] for x in patterns)==list(range(8))
    for row in patterns:
        c=row['corner_code'];assert row['corners']==list(code_bits(c))
        assert row['active_mask_support']==[0,int(bool(c&1)),int(bool(c&2)),int(bool(c))]
    result={'status':'PASS','realizations':len(seen),'family_realizations':dict(sorted(counts.items())),
        'all_exported_networks_checked':True,'all_optima_and_witnesses_checked':True,
        'all_support_patterns_checked':32,'exact_probability_mass_per_family':'1',
        'all_AND_boundary':'Every enumerated n=2 matching has positive physical cost and no marked singleton cycle.',
        'input_sha256':{p.name:sha(p) for p in sorted(directory.glob('*'))},
        'gate_sha256':sha(ROOT/'qa/independent_exact_checks.json'),
        'verification_source_sha256':{p.name:sha(p) for p in [Path(__file__),
            Path(__file__).with_name('independent_exact_checks.py'),
            Path(__file__).with_name('independent_network_oracle.py')]}}
    representation=json.loads((ROOT/'REPRESENTATION.json').read_text())
    assert result['gate_sha256']==representation['current_gate_sha256']
    (ROOT/'qa/independent_exact_exports.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,sort_keys=True))


if __name__=='__main__':main()
