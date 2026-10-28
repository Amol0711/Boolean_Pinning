"""Small-state oracle independent of the production solver.

Physical Boolean evaluation, compatible full-state enumeration, all pin subsets,
nominal/physical support, and elementary marked-slot cycle enumeration. Only
standard-library code is used. No production solver module is imported.
"""
from itertools import combinations, product
from collections import deque


def evaluate(parents, tables, state):
    result = 0
    for i, slots in enumerate(parents):
        corner = sum(((state >> p) & 1) << j for j, p in enumerate(slots))
        result |= int(tables[i][corner]) << i
    return result


def target_integer(target):
    return sum(int(bit) << i for i, bit in enumerate(target))


def success(raw_successors, target, pinned):
    """Reverse reachability from target on every compatible physical state."""
    mask = sum(1 << i for i in pinned)
    compatible = [s for s in range(len(raw_successors)) if (s & mask) == (target & mask)]
    reverse = {s: [] for s in compatible}
    for s in compatible:
        out = (raw_successors[s] & ~mask) | (target & mask)
        reverse[out].append(s)
    reached = {target}
    queue = deque([target])
    while queue:
        state = queue.popleft()
        for predecessor in reverse[state]:
            if predecessor not in reached:
                reached.add(predecessor)
                queue.append(predecessor)
    return len(reached) == len(compatible), len(compatible)


def all_pin_sets(parents, tables, target, baseline=()):
    n = len(parents)
    t = target_integer(target)
    raw = [evaluate(parents, tables, s) for s in range(1 << n)]
    if raw[t] != t:
        raise ValueError('Target is not fixed')
    base = tuple(sorted(baseline))
    free = [i for i in range(n) if i not in base]
    records = []
    optimum = None
    transitions = 0
    for cardinality in range(len(free) + 1):
        for additional in combinations(free, cardinality):
            passed, states = success(raw, t, base + additional)
            transitions += states
            records.append({'additional': additional, 'pins': tuple(sorted(base + additional)),
                            'success': passed, 'states': states})
            if passed and optimum is None:
                optimum = cardinality
    assert optimum is not None
    return {'optimum': optimum, 'raw_successors': raw, 'records': records,
            'compatible_state_checks': transitions}


def support(parents, tables, target, baseline=(), mode='slot'):
    """Brute force over local nominal or consistent physical corners."""
    n = len(parents)
    active = set(range(n)) - set(baseline)
    while True:
        new = set()
        for i in active:
            slots = parents[i]
            if mode == 'slot':
                free = [j for j,p in enumerate(slots) if p in active]
                corner0 = sum(target[p] << j for j,p in enumerate(slots))
                for bits in product((0,1), repeat=len(free)):
                    corner = corner0
                    for j, value in zip(free,bits):
                        corner = (corner & ~(1<<j)) | (value<<j)
                    if tables[i][corner] != target[i]:
                        new.add(i); break
            elif mode == 'physical':
                free = sorted(set(slots) & active)
                for bits in product((0,1),repeat=len(free)):
                    values = dict(zip(free,bits))
                    corner=sum(values.get(p,target[p])<<j for j,p in enumerate(slots))
                    if tables[i][corner] != target[i]:
                        new.add(i);break
            else:
                raise ValueError(mode)
        if new == active:
            return sorted(active)
        assert new <= active
        active = new


def marker(parents, tables, target, receiver, slot):
    corner=sum(target[p]<<j for j,p in enumerate(parents[receiver]))
    return tables[receiver][corner ^ (1<<slot)] != target[receiver]


def marked_slot_cycles(parents,tables,target,baseline=()):
    """Count vertex-simple directed cycles, distinguishing nominal parallel slots.

    Start at the smallest vertex on a cycle. For each vertex sequence multiply
    the number of marked slots for every chosen edge. A repeated-parent OR
    self-loop consequently contributes two cycles, not one.
    """
    nodes=set(range(len(parents)))-set(baseline)
    edge={i:{} for i in nodes}
    for receiver in nodes:
        for j,parent in enumerate(parents[receiver]):
            if parent in nodes and marker(parents,tables,target,receiver,j):
                edge[parent][receiver]=edge[parent].get(receiver,0)+1
    total=0
    for root in sorted(nodes):
        def walk(v, visited, multiplicity):
            subtotal=0
            for w,number in edge[v].items():
                if w==root:
                    subtotal+=multiplicity*number
                elif w>root and w not in visited:
                    subtotal+=walk(w,visited|{w},multiplicity*number)
            return subtotal
        total+=walk(root,{root},1)
    return total


def graph_sccs(nodes,edges):
    """SCCs from independent Floyd-Warshall mutual reachability."""
    nodes=sorted(nodes);reach={v:{v} for v in nodes}
    for a,b in edges:reach[a].add(b)
    for pivot in nodes:
        for a in nodes:
            if pivot in reach[a]:reach[a] |= reach[pivot]
    unseen=set(nodes);comps=[]
    while unseen:
        root=min(unseen)
        C={v for v in unseen if v in reach[root] and root in reach[v]}
        comps.append(sorted(C));unseen-=C
    return comps


def certificates(parents,tables,target,baseline=(),mode='slot'):
    R=support(parents,tables,target,baseline,mode)
    Rset=set(R);edges=[(p,i) for i in R for p in parents[i] if p in Rset]
    G=all(sum(p in Rset for p in parents[i])==1 for i in R)
    Gp=True;cyclic=0
    for C in graph_sccs(R,edges):
        Cset=set(C)
        iscycle=len(C)>1 or any(a==b==C[0] for a,b in edges)
        if not iscycle:continue
        cyclic+=1
        for i in C:
            slots=[j for j,p in enumerate(parents[i]) if p in Cset]
            if len(slots)!=1 or not marker(parents,tables,target,i,slots[0]):
                Gp=False
    return {'residual':R,'G':G,'G_prime':Gp,'cyclic_components':cyclic,
            'marked_cycles':marked_slot_cycles(parents,tables,target,tuple(sorted(set(baseline)|(set(range(len(parents)))-Rset)))),
            'global_marked_cycles':marked_slot_cycles(parents,tables,target,baseline)}
