"""Independent finite checks for the JointRules planted joint-corner law extension.

Run from the archive root: python verification/independent_exact_checks.py
Reference answers use labeled matching enumeration, integer law weights, direct
support evaluation and an independently implemented physical-state oracle.
No production validation or solver modules are imported.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import itertools as it
import json
import math
from pathlib import Path
import random

import numpy as np
from independent_network_oracle import (exact_cost, export_rows,
    graph_certificate, baseline_refined_certificate, literal_support,
    physical_successors, succeeds)

ROOT = Path(__file__).resolve().parents[1]


def multiply(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def code_bits(code):
    return tuple((code >> j) & 1 for j in range(3))


def support_image(parents, codes, mask):
    answer = 0
    for v, ((p, q), code) in enumerate(zip(parents, codes)):
        a, b = (mask >> p) & 1, (mask >> q) & 1
        # Independently free nominal slots, even when p == q.
        if (a and code & 1) or (b and code & 2) or (a and b and code & 4):
            answer |= 1 << v
    return answer


def cycles(parents, codes):
    """Actual marked slot cycles, canonically rotated by smallest vertex."""
    n = len(parents)
    result = []
    for length in range(1, n+1):
        for vertices in it.permutations(range(n), length):
            if vertices[0] != min(vertices):
                continue
            choices = []
            for j, receiver in enumerate(vertices):
                source = vertices[(j-1) % length]
                choices.append([s for s in (0, 1)
                                if parents[receiver][s] == source and
                                ((codes[receiver] >> s) & 1)])
            for slots in it.product(*choices):
                result.append((vertices, slots))
    return result


def rational(value):
    value = Fraction(value)
    return {'numerator': value.numerator, 'denominator': value.denominator}


def fixture_check(name, degrees, laws):
    n = len(degrees)
    assert sum(degrees) == 2*n and len(laws) == n
    denominators = [sum(law.values()) for law in laws]
    law_denominator = math.prod(denominators)
    # Enumerate labeled stubs explicitly; retain equal physical-parent copies.
    stubs = tuple(v for v, d in enumerate(degrees) for _ in range(d))
    matchings = Counter(tuple(stubs[j] for j in perm)
                        for perm in it.permutations(range(2*n)))
    matching_count = math.factorial(2*n)
    assert sum(matchings.values()) == matching_count
    denominator = law_denominator * matching_count
    support_totals = [0] * (1 << n)
    incidence_totals = [[0]*n for _ in range(1 << n)]
    cycle_totals = [0] * (n+1)
    factorial_totals = Counter()
    physical_realizations = repeated_parent_realizations = 0
    total_weight = 0
    law_products = list(it.product(*[list(sorted(law.items())) for law in laws]))
    for flat, matching_weight in sorted(matchings.items()):
        parents = tuple(zip(flat[::2], flat[1::2]))
        for selected in law_products:
            codes = tuple(x[0] for x in selected)
            weight = matching_weight * math.prod(x[1] for x in selected)
            total_weight += weight
            physical_realizations += 1
            repeated_parent_realizations += any(a == b for a, b in parents)
            for mask in range(1 << n):
                if support_image(parents, codes, mask) != mask:
                    continue
                support_totals[mask] += weight
                for v, (a, b) in enumerate(parents):
                    if mask >> v & 1 and mask >> a & 1 and mask >> b & 1:
                        incidence_totals[mask][v] += weight
            actual_cycles = cycles(parents, codes)
            for vertices, _ in actual_cycles:
                cycle_totals[len(vertices)] += weight
            loops = [c for c in actual_cycles if len(c[0]) == 1]
            factorial_totals['loop_factorial_2'] += len(loops)*(len(loops)-1)*weight
            for a, b in it.permutations(actual_cycles, 2):
                if not (set(a[0]) & set(b[0])):
                    factorial_totals[f'disjoint_{len(a[0])}_{len(b[0])}'] += weight
    assert total_weight == denominator
    singleton_weights = [sum(w*((c & 1)+((c >> 1) & 1)) for c,w in law.items())
                         for law in laws]
    nonzero_weights = [sum(w for c,w in law.items() if c) for law in laws]
    q12 = [Fraction(sum(w for c,w in law.items() if c & 3 == 3), den)
           for law,den in zip(laws,denominators)]
    s = [Fraction(w,d) for w,d in zip(singleton_weights,denominators)]
    r = [Fraction(w,d) for w,d in zip(nonzero_weights,denominators)]
    a_factors = [[0,sw,rw] for sw,rw in zip(singleton_weights,nonzero_weights)]
    b_factors = [[d,2*d-sw,d-rw] for d,sw,rw in
                 zip(denominators,singleton_weights,nonzero_weights)]
    fixed_checks = incidence_checks = 0
    for mask in range(1 << n):
        t = sum(d for v,d in enumerate(degrees) if mask >> v & 1)
        factors = [a_factors[v] if mask >> v & 1 else b_factors[v] for v in range(n)]
        poly = [1]
        for f in factors:
            poly = multiply(poly,f)
        assert support_totals[mask]*math.comb(2*n,t) == poly[t]*matching_count
        fixed_checks += 1
        for v in range(n):
            if not mask >> v & 1:
                assert incidence_totals[mask][v] == 0
                continue
            poly = [1]
            for j,f in enumerate(factors):
                poly = multiply(poly,[0,0,nonzero_weights[j]] if j == v else f)
            assert incidence_totals[mask][v]*math.comb(2*n,t) == poly[t]*matching_count
            incidence_checks += 1
    cycle_checks = 0
    weighted_singletons = [d*x for d,x in zip(degrees,s)]
    def distinct_sum(length):
        return sum((math.prod(weighted_singletons[j] for j in seq)
                    for seq in it.permutations(range(n),length)), Fraction(0))
    for length in range(1,n+1):
        expected = distinct_sum(length)/(length*math.perm(2*n,length))
        assert Fraction(cycle_totals[length],denominator) == expected
        cycle_checks += 1
    for l1 in range(1,n+1):
        for l2 in range(1,n+1-l1):
            expected = distinct_sum(l1+l2)/(l1*l2*math.perm(2*n,l1+l2))
            assert Fraction(factorial_totals[f'disjoint_{l1}_{l2}'],denominator) == expected
            cycle_checks += 1
    expected_loop_factorial = (distinct_sum(2) +
        2*sum(d*(d-1)*p for d,p in zip(degrees,q12)))/math.perm(2*n,2)
    assert Fraction(factorial_totals['loop_factorial_2'],denominator) == expected_loop_factorial
    cycle_checks += 1
    lam = sum(weighted_singletons)/(2*n)
    rho = sum(d*x for d,x in zip(degrees,r))/(2*n)
    incidence_bound_checks = []
    admissible = max(lam,rho) < 1
    if admissible:
        b = (1+max(lam,rho))/2
        D = max(degrees)
        h = len(set((d,tuple(sorted(law.items()))) for d,law in zip(degrees,laws)))
        gamma = (math.log(1/float(b))-1+float(b))/D
        for k in range(1,n+1):
            actual = Fraction(sum(sum(incidence_totals[m]) for m in range(1 << n)
                                  if m.bit_count() == k), denominator)
            Lk = 3*math.sqrt(D*k+1)*math.comb(k+h-1,h-1)*math.exp(-gamma*k)
            upper = float(max(r))*D*D*k*k*Lk/(4*float(b)*n)
            assert float(actual) <= upper + 1e-12
            incidence_bound_checks.append({'k':k,'exact_expectation':rational(actual),'bound':upper})
    return {'name':name,'n':n,'degrees':degrees,'law_weights':laws,
            'labeled_matchings':matching_count,'distinct_parent_assignments':len(matchings),
            'physical_rule_realizations':physical_realizations,
            'repeated_parent_realizations':repeated_parent_realizations,
            'fixed_set_coefficient_checks':fixed_checks,
            'distinguished_row_incidence_checks':incidence_checks,
            'cycle_and_factorial_moment_checks':cycle_checks,
            'lambda':rational(lam),'rho':rational(rho),'strict_condition':admissible,
            'zero_safe_incidence_bound_checks':incidence_bound_checks}


def physical_checks():
    counts = Counter()
    def check(parents,codes, exhaustive_pins):
        n = len(parents)
        corners = np.array([code_bits(c) for c in codes],dtype=np.uint8)
        parents = np.array(parents,dtype=np.int64)
        residual,_ = literal_support(parents,corners)
        rows = export_rows(parents,corners,list(range(n)))
        residual_rows = export_rows(parents,corners,residual)
        graph = graph_certificate(residual_rows)
        global_result = exact_cost(rows,list(range(n)))
        residual_result = exact_cost(residual_rows,residual)
        assert global_result['complete'] and residual_result['complete']
        optimum = global_result['cost']
        assert optimum == residual_result['cost'] <= len(residual)
        local_sum = sum(exact_cost(residual_rows,c['nodes'])['cost'] for c in graph['components'])
        assert local_sum == optimum
        marked = len(cycles(parents.tolist(),codes))
        if graph['G_prime']:
            assert optimum == marked == graph['cyclic_component_count']
            counts['G_prime_cost_equalities'] += 1
        assert not graph['G'] or graph['G_prime']
        certificate = baseline_refined_certificate(residual_rows)
        assert certificate['lower'] <= optimum <= certificate['upper'] <= len(residual)
        successor = physical_successors(rows,list(range(n)))
        assert succeeds(successor,sum(1 << v for v in certificate['pins']))
        if exhaustive_pins:
            residual_successor = physical_successors(residual_rows,residual)
            for mask in range(1 << n):
                reduced_mask = sum(1 << j for j,v in enumerate(residual) if mask >> v & 1)
                assert succeeds(successor,mask) == succeeds(residual_successor,reduced_mask)
                counts['pin_set_capture_equalities'] += 1
                # Independent literal trajectory check of the state-graph oracle.
                literal_success = True
                for start in range(1 << n):
                    if start & mask:
                        continue
                    state,seen = start,set()
                    while state and state not in seen:
                        seen.add(state)
                        image = 0
                        for v,(p,q) in enumerate(parents):
                            index = ((state >> int(p)) & 1)+2*((state >> int(q)) & 1)
                            if index and (codes[v] >> (index-1)) & 1:
                                image |= 1 << v
                        state = image & ~mask
                    literal_success &= state == 0
                    counts['literal_compatible_trajectories'] += 1
                assert succeeds(successor,mask) == literal_success
        counts['networks'] += 1
        counts['global_optima'] += 1
        counts['component_optima'] += len(graph['components'])
    for flat in it.product(range(2),repeat=4):
        for codes in it.product(range(8),repeat=2):
            check(tuple(zip(flat[::2],flat[1::2])),codes,True)
    counts['exhaustive_two_node_networks'] = counts['networks']
    rng = random.Random(2026101504)
    for _ in range(300):
        n = rng.choice([3,4,5,6])
        check([(rng.randrange(n),rng.randrange(n)) for _ in range(n)],
              [rng.randrange(8) for _ in range(n)],False)
    counts['additional_seeded_networks'] = 300
    return dict(counts)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'qa/independent_exact_checks.json')
    args = parser.parse_args()
    fixtures = [
      ('one_vertex_dependent',[2],[{0:2,3:1,4:1}]),
      ('asymmetric_degree_correlated',[1,3],[{0:2,5:3,6:1,4:2},{0:3,7:1,3:2,4:2}]),
      ('zero_singleton_type',[1,3],[{0:1,5:1},{0:2,4:3}]),
      ('zero_degree_three_types',[0,2,4],[{0:4,5:3,6:1},{0:5,7:1,3:2},{0:1,4:1}]),
      ('all_eight_dependent_corners',[2,2,2],[dict(enumerate([8,1,2,3,4,5,6,7]))]*3),
      ('lambda_zero_admissible',[2,2,2],[{0:1,4:1}]*3),
      ('constant_zero',[1,2,3],[{0:1}]*3),
      ('all_AND_negative_condition',[2,2,2],[{4:1}]*3),
      ('independent_Bernoulli_quarter',[1,3],[{c:3**(3-c.bit_count()) for c in range(8)}]*2),
    ]
    protocol_path = ROOT/'validation/Protocol_JointRules.json'
    protocol = json.loads(protocol_path.read_text())
    for family,spec in sorted(protocol['families'].items()):
        degrees,laws = [],[]
        for row in spec['types']:
            count = 2*Fraction(*row['fraction'])
            assert count.denominator == 1
            degrees.extend([row['degree']]*int(count))
            laws.extend([{c:w for c,w in enumerate(row['weights']) if w}]*int(count))
        fixtures.append((family,degrees,laws))
    records = [fixture_check(*fixture) for fixture in fixtures]
    assert not next(x for x in records if x['name']=='all_AND_negative_condition')['strict_condition']
    result = {'status':'PASS','scope':'Exact finite-state and rational-arithmetic validation checks.',
              'fixtures':records,'physical_checks':physical_checks(),
              'protocol_sha256':hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
              'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                  [Path(__file__),Path(__file__).with_name('independent_network_oracle.py')]}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'fixtures':len(records),
                      'physical_checks':result['physical_checks']},sort_keys=True))


if __name__ == '__main__':
    main()
