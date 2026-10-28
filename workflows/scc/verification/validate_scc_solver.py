"""Deterministic independent small-network validation.

Run from the package root:
  python verification/validate_scc_solver.py --output results/solver_validation.json
Only the production constructor and public solver are imported. All state-map,
feasibility, exact-optimum, support and cycle oracles are independently written.
No NumPy, external validation package, model files, or network access is required.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import random
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from scc_control import BooleanNetwork,solve_network,evaluate_state
from independent_oracle import all_pin_sets,evaluate,target_integer,success,certificates


def fixtures():
    data=[
      {'name':'repeated_self_xor_physical_zero','parents':[[0,0]],'tables':[[0,1,1,0]],'target':[0],'baseline':[], 'expected_cost':0,'expected_slot_residual':[0],'expected_physical_residual':[],'expected_G':False,'expected_G_prime':False,'expected_marked_cycles':2},
      {'name':'repeated_self_or_parallel_marked_slots','parents':[[0,0]],'tables':[[0,1,1,1]],'target':[0],'baseline':[],'expected_cost':1,'expected_slot_residual':[0],'expected_physical_residual':[0],'expected_G':False,'expected_G_prime':False,'expected_marked_cycles':2},
      {'name':'repeated_self_and_synergy_unmarked','parents':[[0,0]],'tables':[[0,0,0,1]],'target':[0],'baseline':[],'expected_cost':1,'expected_slot_residual':[0],'expected_physical_residual':[0],'expected_G':False,'expected_G_prime':False,'expected_marked_cycles':0},
      {'name':'downstream_nonlinear_Gfalse_Gprimetrue','parents':[[0,2],[0,0],[]],'tables':[[0,1,0,1],[0,1,1,1],[0]],'target':[0,0,0],'baseline':[],'expected_cost':1,'expected_slot_residual':[0,1],'expected_physical_residual':[0,1],'expected_G':False,'expected_G_prime':True,'expected_marked_cycles':1},
      {'name':'ignored_self_slot_then_frozen','parents':[[0,1],[]],'tables':[[0,0,1,1],[0]],'target':[0,0],'baseline':[],'expected_cost':0,'expected_slot_residual':[],'expected_physical_residual':[],'expected_G':True,'expected_G_prime':True,'expected_marked_cycles':0},
      {'name':'two_node_complex_cost_two_budget_case','parents':[[0,1],[0,1]],'tables':[[0,1,1,1],[0,1,1,1]],'target':[0,0],'baseline':[],'expected_cost':2,'expected_slot_residual':[0,1],'expected_physical_residual':[0,1],'expected_G':False,'expected_G_prime':False,'expected_marked_cycles':3},
      {'name':'mixed_target_baseline_and_constant','parents':[[],[2],[2]],'tables':[[1],[1,0],[0,1]],'target':[1,0,1],'baseline':[2],'expected_cost':0,'expected_slot_residual':[],'expected_physical_residual':[],'expected_G':True,'expected_G_prime':True,'expected_marked_cycles':0},
    ]
    # Eight nominal slots, each of four physical coordinates appearing twice.
    parents=[[0,1,2,3,0,1,2,3],[0,0,1,2,3],[1,3],[]]
    target=[1,0,1,1]
    tables=[]
    for i,slots in enumerate(parents):
        tab=[((mask.bit_count()+i)%2) for mask in range(1<<len(slots))]
        fixed=sum(target[p]<<j for j,p in enumerate(slots));tab[fixed]=target[i]
        tables.append(tab)
    data.append({'name':'high_arity_eight_slots_mixed_target_baseline','parents':parents,'tables':tables,'target':target,'baseline':[3]})
    return data


def remaining_optimum(successful,pins):
    mask=sum(1<<i for i in pins)
    return min((m&~mask).bit_count() for m in successful if (m&mask)==mask)


def check_result(result,case,baseline,oracle,expected,mode,stats):
    p,t,x=case['parents'],case['tables'],case['target']
    independent=certificates(p,t,x,baseline,mode)
    assert result['residual']==independent['residual'],('support',case['name'],baseline,mode,result,independent)
    assert result['G']==independent['G'],('G',case['name'],baseline)
    assert result['G_prime']==independent['G_prime'],('Gprime',case['name'],baseline)
    assert result['cyclic_component_count']==independent['cyclic_components']
    assert not result['G'] or result['G_prime']
    assert result['lower']<=expected<=result['upper'],('bounds',case['name'],baseline,result,expected)
    if result['exact_cost'] is not None:assert result['exact_cost']==expected
    pins=tuple(result['feasible_pins'])
    assert len(set(pins))==len(pins)==result['upper']
    assert not set(pins)&set(baseline)
    assert success(oracle['raw_successors'],target_integer(x),tuple(baseline)+pins)[0],('witness',case['name'],baseline,pins)
    if mode=='slot':assert independent['marked_cycles']==independent['global_marked_cycles']
    if result['G_prime']:
        assert expected==independent['marked_cycles']==result['certificate_cycle_count'],('cyclecost',case['name'],baseline,independent,result)
        stats['G_prime_cost_cycle_equalities']+=1
    else:assert result['certificate_cycle_count'] is None
    stats['solver_calls']+=1
    stats['bound_checks']+=1
    stats['feasible_witness_checks']+=1
    stats['support_and_certificate_checks']+=1
    return independent


def check_case(case,stats,all_baselines=True,modes=('slot',)):
    p,t,x=case['parents'],case['tables'],case['target'];base=tuple(case.get('baseline',[]))
    net=BooleanNetwork(p,t,x)
    oracle=all_pin_sets(p,t,x,base)
    stats['networks']+=1
    stats['full_pin_subset_oracle_checks']+=len(oracle['records'])
    stats['compatible_state_successor_checks']+=oracle['compatible_state_checks']
    stats['raw_physical_state_checks']+=1<<net.n
    # Check public physical evaluator independently, including time-zero clamping.
    for state in range(1<<net.n):
        assert evaluate_state(net,state)==oracle['raw_successors'][state]
        mask=sum(1<<i for i in base);target=target_integer(x)
        clamped=(state&~mask)|(target&mask)
        expected=(evaluate(p,t,clamped)&~mask)|(target&mask)
        assert evaluate_state(net,state,base)==expected
    successful=[sum(1<<i for i in r['pins']) for r in oracle['records'] if r['success']]
    base_sets=[r['pins'] for r in oracle['records']] if all_baselines else [base]
    selected=[]
    for mode in modes:
      for K in base_sets:
        expected=remaining_optimum(successful,K)
        result=solve_network(net,baseline=K,support_mode=mode,max_component_size=12,max_pins=12,max_work=10_000_000)
        assert result['exact_cost'] is not None,('unexpected incomplete',case['name'],K,result)
        independent=check_result(result,case,K,oracle,expected,mode,stats)
        if K==base:
            selected.append({'mode':mode,'optimum':expected,'residual':independent['residual'],'G':independent['G'],'G_prime':independent['G_prime'],'marked_cycles':independent['marked_cycles']})
            if mode=='slot':
                for key,value in [('expected_cost',expected),('expected_slot_residual',independent['residual']),('expected_G',independent['G']),('expected_G_prime',independent['G_prime']),('expected_marked_cycles',independent['marked_cycles'])]:
                    if key in case:assert case[key]==value,(case['name'],key,value)
            elif 'expected_physical_residual' in case:
                assert case['expected_physical_residual']==independent['residual']
    return selected,oracle


def empty_stats():
    return {k:0 for k in ['networks','solver_calls','full_pin_subset_oracle_checks','compatible_state_successor_checks','raw_physical_state_checks','bound_checks','feasible_witness_checks','support_and_certificate_checks','G_prime_cost_cycle_equalities']}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE.parent/'results/solver_validation.json')
    args=parser.parse_args()
    out={'schema_version':'1.0','status':'PASS','scope':'Finite explicitly enumerated validation suite, not an exhaustive verification of arbitrary network sizes.','oracle':'Independent physical full-state maps and every additional pin subset; independent support/SCC/marked-slot-cycle computations.','random_seed':2026101501,'suites':{}}
    stats=empty_stats();seen=0
    for slots in itertools.product(range(2),repeat=4):
      p=[list(slots[:2]),list(slots[2:])]
      for a,b in itertools.product(range(8),repeat=2):
        t=[[0,(code>>0)&1,(code>>1)&1,(code>>2)&1] for code in (a,b)]
        case={'name':f'n2_{seen}','parents':p,'tables':t,'target':[0,0]}
        check_case(case,stats);seen+=1
    assert seen==1024
    out['suites']['exhaustive_two_node_two_slot']=stats
    print('PASS 1024 planted two-node networks, all pin subsets',flush=True)
    stats=empty_stats();rng=random.Random(out['random_seed']);digest=hashlib.sha256();by_n={}
    for i in range(300):
        n=3+i%6;target=[rng.randrange(2) for _ in range(n)];p=[];t=[]
        for v in range(n):
            arity=rng.randrange(5);row=[rng.randrange(n) for _ in range(arity)]
            tab=[rng.randrange(2) for _ in range(1<<arity)]
            tab[sum(target[parent]<<j for j,parent in enumerate(row))]=target[v]
            p.append(row);t.append(tab)
        baseline=sorted(rng.sample(range(n),rng.randrange(min(2,n)+1)))
        case={'name':f'random_{i}','parents':p,'tables':t,'target':target,'baseline':baseline}
        digest.update(json.dumps(case,sort_keys=True,separators=(',',':')).encode())
        check_case(case,stats,modes=('slot','physical'))
        by_n[str(n)]=by_n.get(str(n),0)+1
    stats['by_network_size']=by_n;stats['case_stream_sha256']=digest.hexdigest()
    out['suites']['random_arbitrary_rules_and_baselines']=stats
    print('PASS 300 arbitrary-rule networks n=3..8, all pin subsets, both support modes',flush=True)
    fixed=fixtures();(HERE/'fixtures.json').write_text(json.dumps(fixed,indent=2,sort_keys=True)+'\n')
    stats=empty_stats();records=[]
    for case in fixed:
        results,oracle=check_case(case,stats,modes=('slot','physical'))
        records.append({'name':case['name'],'results':results})
    out['suites']['boundary_fixtures']={'coverage':stats,'cases':records}
    print('PASS eight boundary fixtures, including high arity and mixed target',flush=True)
    # Distinct truncation mechanisms, including interruption during candidates.
    case=next(c for c in fixed if c['name']=='two_node_complex_cost_two_budget_case')
    net=BooleanNetwork(case['parents'],case['tables'],case['target']);oracle=all_pin_sets(case['parents'],case['tables'],case['target'])
    budget_specs=[('component_size',{'max_component_size':1},'component_size_cap'),('pin_count',{'max_pins':0},'pin_count_cap'),('work_preflight',{'max_work':0},'work_cap_successor_construction'),('work_search',{'max_work':8},'work_cap_candidate_search')]
    stats=empty_stats();records=[]
    for label,options,reason in budget_specs:
        result=solve_network(net,**options)
        check_result(result,case,(),oracle,2,'slot',stats)
        reasons=[c['reason'] for c in result['components']]
        assert reason in reasons,(label,reasons,reason)
        assert result['exact_cost'] is None,(label,result)
        records.append({'case':label,'options':options,'lower':result['lower'],'upper':result['upper'],'reasons':reasons,'feasible_pins':result['feasible_pins']})
    out['suites']['incomplete_budgets']={'coverage':stats,'cases':records}
    print('PASS size, pin-count and work-budget truncation bounds/witnesses',flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print('PASS deterministic validation record written',flush=True)

if __name__=='__main__':main()
