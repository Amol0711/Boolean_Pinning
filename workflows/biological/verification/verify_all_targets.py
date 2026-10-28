"""Independent all-target biological verification (Python standard library).

  python verification/verify_all_targets.py --production results/all_targets.json

Own BNET grammar/local tables (independent_core), arbitrary-precision bitset
fixed-point enumeration, exact essential-dependency graph, Floyd SCCs, reverse
reachability for every tested local pin set, and all-compatible-residual-state
checks of the reported global witnesses. No production solver is imported.
"""
from pathlib import Path
from itertools import product,combinations
from collections import Counter
import argparse,hashlib,json,sys
import independent_core as IC

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read_asts(path):
    out={}
    for line in Path(path).read_text().splitlines():
        line=line.strip()
        if not line or line.startswith('#') or line.startswith('targets,'):continue
        name,text=line.split(',',1);out[name.strip()]=IC.parse_expression(text)
    return out


def force_values(rows,inputs,assignment):
    fixed={n:assignment[n] for n in inputs}
    while True:
        changed=False
        for name,scope,table in rows:
            if name in fixed:continue
            possible={v for bits,v in table.items() if all(p not in fixed or fixed[p]==b for p,b in zip(scope,bits))}
            if len(possible)==1:fixed[name]=possible.pop();changed=True
        if not changed:return fixed


def bit_pattern(j,count):
    width=1<<(j+1);out=((1<<(1<<j))-1)<<(1<<j)
    while width<count:out|=out<<width;width*=2
    return out


def bit_eval(tree,values,mask):
    op=tree[0]
    if op=='variable':return values[tree[1]]
    if op=='constant':return mask if tree[1] else 0
    if op=='not':return mask^bit_eval(tree[1],values,mask)
    a=bit_eval(tree[1],values,mask);b=bit_eval(tree[2],values,mask)
    return a&b if op=='and' else a|b


def enumerate_fixed_points(rows,regulated,inputs,asts,assignment,order):
    forced=force_values(rows,inputs,assignment)
    free=sorted(set(regulated)-set(forced))
    N=1<<len(free);mask=(1<<N)-1
    values={n:mask if v else 0 for n,v in forced.items()}
    values.update({n:bit_pattern(j,N) for j,n in enumerate(free)})
    good=mask
    # Every regulated equation is tested, including all previously forced rows.
    for name in regulated:good&=mask^(bit_eval(asts[name],values,mask)^values[name])
    targets=[]
    while good:
        bit=good&-good;state=bit.bit_length()-1;good^=bit
        target=forced.copy();target.update({n:(state>>j)&1 for j,n in enumerate(free)})
        assert IC.step(rows,target)==target
        targets.append(target)
    targets.sort(key=lambda target:tuple(target[n] for n in order))
    return forced,free,N,targets


def functional_rows(rows,target,residual):
    R=set(residual);out=[]
    for name,scope,table in rows:
        if name not in R:continue
        active=sorted(set(scope)&R)
        restricted=[]
        for state in range(1<<len(active)):
            corner=target.copy();corner.update({n:target[n]^((state>>j)&1) for j,n in enumerate(active)})
            restricted.append(table[tuple(corner[p] for p in scope)]^target[name])
        essential=[j for j in range(len(active)) if any(restricted[s]!=restricted[s^(1<<j)] for s in range(len(restricted)))]
        parents=[active[j] for j in essential]
        compressed=[restricted[sum(((s>>j)&1)<<old for j,old in enumerate(essential))] for s in range(1<<len(essential))]
        out.append({'name':name,'parents':parents,'table':compressed})
    return sorted(out,key=lambda r:r['name'])


def sccs(nodes,rows):
    nodes=sorted(nodes);reach={v:{v} for v in nodes};edges=[]
    for row in rows:
        for p in row['parents']:reach[p].add(row['name']);edges.append((p,row['name']))
    for k in nodes:
        for v in nodes:
            if k in reach[v]:reach[v]|=reach[k]
    rem=set(nodes);components=[]
    while rem:
        v=min(rem);C=sorted(w for w in rem if w in reach[v] and v in reach[w]);rem-=set(C)
        components.append({'names':C,'cyclic':len(C)>1 or (v,v) in edges})
    return components


def minimum_local(rows,target,nodes,eligible=None):
    nodes=sorted(nodes);eligible=nodes if eligible is None else sorted(set(eligible)&set(nodes))
    allowed=[nodes.index(n) for n in eligible]
    successors=IC.residual_successors(rows,target,nodes)
    checked=0;statechecks=0;levels=[]
    for size in range(len(allowed)+1):
        minimizers=[];tested=0
        for K in combinations(allowed,size):
            r=IC.controlled_graph(successors,K);tested+=1;checked+=1;statechecks+=r['compatible_states']
            if r['success']:minimizers.append([nodes[j] for j in K])
        levels.append({'cardinality':size,'tested':tested,'successful':len(minimizers)})
        if minimizers:return {'cost':size,'minimizers':minimizers,'levels':levels,'candidate_count':checked,'compatible_states_checked':statechecks}
    return {'cost':None,'minimizers':[],'levels':levels,'candidate_count':checked,'compatible_states_checked':statechecks}


def global_witness(rows,target,residual,pins):
    residual=set(residual);pins=set(pins);free=sorted(residual-pins)
    # Verify all unpinned rows outside residual remain at target on the whole cube.
    for name,scope,table in rows:
        if name in residual:continue
        for bits,value in table.items():
            if all(p in residual or bit==target[p] for p,bit in zip(scope,bits)):
                assert value==target[name],('noninvariant',name)
    successors=IC.residual_successors(rows,target,free)
    answer=IC.controlled_graph(successors,())
    assert answer['success']
    return {'compatible_residual_states':answer['compatible_states'],'maximum_residual_hitting_time':answer['maximum_hitting_time'],'success':True}


def discover(rows,assignment,order):
    state={n:assignment.get(n,0) for n in order};seen={};trail=[]
    for step in range(10001):
        key=tuple(state[n] for n in order)
        if key in seen:
            first=seen[key];cycle=trail[first:]
            return state if len(cycle)==1 else None,{'first_repeated_index':first,'repeat_update':step,'cycle_length':len(cycle)}
        seen[key]=step;trail.append(key);state=IC.step(rows,state)
    raise RuntimeError('Orbit budget exceeded')


def expected_records():
    protocol=json.loads((ROOT/'larger_networks/protocol/contexts.json').read_text())
    models={};out=[];totals=Counter()
    for mid in ['cd4_188','macrophage_321']:
        path=ROOT/'larger_networks/data'/mid/'model.bnet';rows,regulated,inputs=IC.read_bnet(path)
        order=sorted(regulated)+sorted(inputs);asts=read_asts(path)
        assert sha(path)==protocol['models'][mid]['bnet_sha256']
        models[mid]={'node_order':order,'source_sha256':sha(path)}
        for ctx in [c for c in protocol['contexts'] if c['model_id']==mid]:
            assignment=ctx['input_assignment'];cid=ctx['id']
            forced,free,N,targets=enumerate_fixed_points(rows,regulated,inputs,asts,assignment,order)
            selected,orbit=discover(rows,assignment,order)
            rec={'id':cid,'model_id':mid,'input_assignment':assignment,'forced_values':forced,'free_names':free,'enumerated_assignments':N,'fixed_point_count':len(targets),'zero_orbit':orbit,'fixed_points':[]}
            for target in targets:
                residual=sorted(IC.support_trace(rows,target,inputs)[-1]);fr=functional_rows(rows,target,residual);components=sccs(residual,fr)
                cost=0;pins=[];forbidden_cost=0;forbidden_possible=True;pieces=[]
                for comp in components:
                    local=minimum_local(rows,target,comp['names'])
                    assert local['cost'] is not None
                    cost+=local['cost'];pins+=local['minimizers'][0]
                    totals['local_candidates']+=local['candidate_count'];totals['local_compatible_states']+=local['compatible_states_checked']
                    forbidden=None
                    if 'v_proliferation' in comp['names']:
                        forbidden=minimum_local(rows,target,comp['names'],set(comp['names'])-{'v_proliferation'})
                        totals['restricted_local_candidates']+=forbidden['candidate_count'];totals['restricted_local_compatible_states']+=forbidden['compatible_states_checked']
                        if forbidden['cost'] is None:forbidden_possible=False
                        else:forbidden_cost+=forbidden['cost']
                    else:forbidden_cost+=local['cost']
                    pieces.append(dict(comp,**local,proliferation_forbidden=forbidden))
                witness=global_witness(rows,target,residual,pins)
                totals['global_witness_compatible_states']+=witness['compatible_residual_states'];totals['targets']+=1
                if mid=='cd4_188' and target['v_proliferation']==0:
                    _,scope,tab=next(r for r in rows if r[0]=='v_proliferation');j=scope.index('v_proliferation')
                    assert all(value==1 for bits,value in tab.items() if bits[j]==1)
                    assert not forbidden_possible
                targetrec={'target': [target[n] for n in order],'target_by_name':target,'selected_zero_orbit':selected==target,'residual_names':residual,'functional_rows':fr,'components':pieces,'exact_cost':cost,'global_pin_names':sorted(pins),'global_witness':witness,'proliferation_forbidden':{'feasible':forbidden_possible,'exact_cost':forbidden_cost if forbidden_possible else None}}
                rec['fixed_points'].append(targetrec)
            if selected is not None:
                saved=json.loads((ROOT/'larger_networks/results'/cid/'certificate.json').read_text())
                match=next(t for t in rec['fixed_points'] if t['selected_zero_orbit'])
                assert match['target_by_name']==saved['target_by_name']
                assert len(match['residual_names'])==saved['residual_size']
                assert match['exact_cost']==saved['additional_cost_bounds'][0]==saved['additional_cost_bounds'][1]
            rec['cost_histogram']={str(k):v for k,v in sorted(Counter(t['exact_cost'] for t in rec['fixed_points']).items())}
            out.append(rec);totals['contexts']+=1;totals['fixedpoint_candidate_assignments']+=N
            print(cid,len(targets),'fixed points; cost histogram',rec['cost_histogram'],flush=True)
    assert [c['fixed_point_count'] for c in out if c['model_id']=='cd4_188']==[6,6,15,6,7,6]
    assert next(c for c in out if c['id']=='MAC_BACTERIA')['fixed_point_count']==0
    return {'models':models,'contexts':out,'counts':dict(totals)}


def compare_production(expected,path):
    production=json.loads(path.read_text())
    assert production['protocol_sha256']==sha(ROOT/'larger_networks/protocol/contexts.json')
    modelrows={}
    for mid,em in expected['models'].items():
        gm=production['models'][mid]
        assert gm['node_order']==em['node_order']
        assert gm['bnet_sha256']==em['source_sha256']
        modelrows[mid]=IC.read_bnet(ROOT/'larger_networks/data'/mid/'model.bnet')[0]
    checks=[];counts=Counter()
    actual={c['id']:c for c in production['contexts']}
    assert set(actual)=={c['id'] for c in expected['contexts']}
    for exp in expected['contexts']:
        got=actual[exp['id']];mid=exp['model_id'];order=expected['models'][mid]['node_order'];rows=modelrows[mid]
        for key in ['fixed_point_count','input_assignment','forced_values','enumerated_assignments','cost_histogram']:assert got[key]==exp[key],(exp['id'],key)
        assert set(got['free_names'])==set(exp['free_names'])
        assert [order[i] for i in got['free_nodes']]==got['free_names']
        assert got['zero_orbit_discovery']['period']==exp['zero_orbit']['cycle_length']
        assert got['zero_orbit_discovery']['cycle_start']==exp['zero_orbit']['first_repeated_index']
        cycle=got['zero_orbit_discovery']['cycle']
        for x,y in zip(cycle,cycle[1:]+cycle[:1]):assert IC.step(rows,dict(zip(order,x)))==dict(zip(order,y))
        def named_rows(recordrows):
            return sorted([{'name':order[r['node']],'parents':[order[i] for i in r['parents']],'table':r['table']} for r in recordrows],key=lambda r:r['name'])
        gtargets={tuple(t['target']):t for t in got['fixed_points']}
        assert len(gtargets)==len(got['fixed_points'])
        assert set(gtargets)=={tuple(t['target']) for t in exp['fixed_points']}
        selected=[]
        for t in exp['fixed_points']:
            g=gtargets[tuple(t['target'])];target=t['target_by_name'];R=t['residual_names']
            assert g['target_by_name']==target and g['selected_zero_orbit']==t['selected_zero_orbit']
            assert g['exact_cost']==t['exact_cost']
            assert set(g['residual_names'])==set(R)
            assert [order[i] for i in g['residual']]==g['residual_names']
            assert named_rows(g['functional_rows'])==t['functional_rows']
            assert all(target['v_'+name]==value for name,value in g['marker_signature'].items())
            if t['selected_zero_orbit']:selected.append(g['target_id'])
            components={tuple(c['names']):c for c in g['components']}
            assert set(components)=={tuple(c['names']) for c in t['components']}
            for ec in t['components']:
                gc=components[tuple(ec['names'])]
                assert [order[i] for i in gc['nodes']]==gc['names']
                assert gc['cyclic']==ec['cyclic'] and gc['exact_cost']==ec['cost']
                assert named_rows(gc['restricted_rows'])==functional_rows(rows,target,ec['names'])
                minimum_names=sorted([sorted(order[i] for i in K) for K in gc['minimum_pinsets']])
                assert minimum_names==sorted(ec['minimizers'])
                assert sorted(gc['minimum_pinsets_by_name'])==minimum_names
                assert gc['all_minimum_pinsets_enumerated'] is True
                assert sorted(gc['pin_names']) in minimum_names
                counts['components_compared']+=1
            pin_names=[order[i] for i in g['global_pins']]
            assert pin_names==g['global_pin_names'] and len(pin_names)==len(set(pin_names))==g['exact_cost']
            assert set(pin_names)<=set(R)
            witness=global_witness(rows,target,R,pin_names)
            assert witness['compatible_residual_states']==g['global_witness_compatible_residual_states']
            assert g['global_witness_exhaustively_verified'] is True
            counts['production_global_witness_states']+=witness['compatible_residual_states']
            eligibility=g['proliferation_forbidden'];ef=t['proliferation_forbidden']
            if mid=='cd4_188':
                assert eligibility['cost']==ef['exact_cost']
                assert (eligibility['status']=='exact')==ef['feasible']
                assert eligibility['target_proliferation']==target['v_proliferation']
                assert eligibility['target_STAT5_HIGH']==target['v_STAT5_HIGH']
                if ef['feasible']:
                    epins=[order[i] for i in eligibility['global_pins']]
                    assert 'v_proliferation' not in epins and len(epins)==ef['exact_cost']
                    ewitness=global_witness(rows,target,R,epins)
                    counts['restricted_production_witness_states']+=ewitness['compatible_residual_states']
                else:
                    assert eligibility['global_pins'] is None
                counts['restricted_eligibility_targets']+=1
            else:assert eligibility['status']=='not_applicable'
            counts['targets_compared']+=1
        assert got['selected_target_id']==(selected[0] if selected else None)
        checks.append(exp['id'])
    return {'production_sha256':sha(path),'contexts_compared':checks,'counts':dict(counts),'checks':['All source hashes and original frozen input assignments','Every fixed-point equation over all remaining assignments','Target ordering, full target maps and raw marker signatures','Exact physical residual and essential functional rows','Full essential-dependency SCC partition and every isolated truth table','All local minimum pin sets, independently rejecting all smaller cardinalities','Every reported full residual global witness from every compatible state','All proliferation-ineligible feasibility/cost values and feasible witnesses','All ten prior selected target vectors, residual sizes and exact costs unchanged','BACTERIA absence of fixed points by all 524288 residual assignments']}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--production',type=Path,default=ROOT/'results/all_targets.json')
    parser.add_argument('--output',type=Path,default=ROOT/'results/independent_biological_validation.json')
    parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args()
    expected=expected_records()
    record={'schema_version':'1.0','status':'PASS','method':'Independent standard-library parser, exhaustive bitset fixed-point equations, essential dependency SCCs, complete tested local pin-cardinality levels and compatible residual-state witness checks.','expected':expected}
    if not args.prepare_only:record['production_comparison']=compare_production(expected,args.production)
    else:record['status']='EXPECTED_DATA_PREPARED_PRODUCTION_COMPARISON_PENDING'
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    print(record['status'],flush=True)

if __name__=='__main__':main()
