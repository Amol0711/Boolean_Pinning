"""Exact local-cube certificates for synchronous, target-pinned Boolean maps.

Only Python's standard library is required. Indices in this API are zero based.
Each truth table uses DISTINCT physical variables; its index is the little-endian
binary encoding of the parent values in the stated order. There is no treatment
of repeated references to one variable as independent inputs.
"""

from dataclasses import dataclass
from itertools import product


@dataclass(frozen=True)
class LocalRule:
    parents: tuple
    table: tuple

    def __post_init__(self):
        if len(set(self.parents)) != len(self.parents):
            raise ValueError("Truth-table parents must be distinct physical variables")
        if any(type(p) is not int or p < 0 for p in self.parents):
            raise ValueError("Parents must be nonnegative integer indices")
        if len(self.table) != 1 << len(self.parents):
            raise ValueError("Truth-table length does not match arity")
        if any(v not in (0, 1) for v in self.table):
            raise ValueError("Truth-table values must be Boolean")

    def evaluate(self, state):
        return self.table[sum(state[p] << j for j, p in enumerate(self.parents))]


@dataclass(frozen=True)
class BooleanModel:
    rules: tuple
    target: tuple

    def __post_init__(self):
        if len(self.rules) != len(self.target) or not self.rules:
            raise ValueError("A nonempty model needs one target bit per row")
        if any(t not in (0, 1) for t in self.target):
            raise ValueError("Target values must be Boolean")
        if any(p >= len(self.rules) for r in self.rules for p in r.parents):
            raise ValueError("Parent index outside the network")
        if tuple(r.evaluate(self.target) for r in self.rules) != self.target:
            raise ValueError("The specified target is not a fixed point")

    @property
    def n(self):
        return len(self.rules)

    def indices(self, nodes):
        result = frozenset(nodes)
        if any(type(p) is not int or not 0 <= p < self.n for p in result):
            raise ValueError("Pin/support index outside the network")
        return result

    def clamp(self, state, pins=()):
        pins = self.indices(pins)
        if len(state) != self.n or any(x not in (0, 1) for x in state):
            raise ValueError("State must be a Boolean vector of network size")
        return tuple(self.target[i] if i in pins else x for i, x in enumerate(state))

    def step(self, state, pins=()):
        """Synchronous update; pins are enforced at time zero and after updates."""
        pins = self.indices(pins)
        state = self.clamp(state, pins)
        return tuple(self.target[i] if i in pins else r.evaluate(state)
                     for i, r in enumerate(self.rules))

    def relabeled(self):
        """Return g(y)=f(y XOR target) XOR target, with the same parent order."""
        transformed = []
        for i, rule in enumerate(self.rules):
            entries = []
            for code in range(1 << len(rule.parents)):
                physical_code = sum((((code >> j) & 1) ^ self.target[p]) << j
                                    for j, p in enumerate(rule.parents))
                entries.append(rule.table[physical_code] ^ self.target[i])
            transformed.append(LocalRule(rule.parents, tuple(entries)))
        return BooleanModel(tuple(transformed), (0,) * self.n)


def local_rule_from_callable(references, function):
    """Compile a rule on nominal references into a physical-variable table.

    `function` receives the nominal input values; repeated references therefore
    receive the same value. Distinct parents retain first-occurrence order.
    """
    parents = tuple(dict.fromkeys(references))
    location = {p: j for j, p in enumerate(parents)}
    entries = []
    for code in range(1 << len(parents)):
        args = [((code >> location[p]) & 1) for p in references]
        value = function(*args)
        if value not in (0, 1):
            raise ValueError("Rule callable returned a non-Boolean value")
        entries.append(int(value))
    return LocalRule(parents, tuple(entries))


def cube_row_table(zero_model, row, support):
    """Exact row restriction to Q_support in target-relative coordinates."""
    support = zero_model.indices(support)
    rule = zero_model.rules[row]
    active = tuple(p for p in rule.parents if p in support)
    location = {p: j for j, p in enumerate(active)}
    entries = []
    for code in range(1 << len(active)):
        full_code = sum((((code >> location[p]) & 1) << j)
                        for j, p in enumerate(rule.parents) if p in support)
        entries.append(rule.table[full_code])
    return active, tuple(entries)


def support_map(zero_model, support, pins=()):
    if any(zero_model.target):
        raise ValueError("Support routines require the relabeled zero-target model")
    pins = zero_model.indices(pins)
    support = zero_model.indices(support)
    return frozenset(i for i in range(zero_model.n) if i not in pins
                     and any(cube_row_table(zero_model, i, support)[1]))


def support_trace(model, pins=()):
    """Descending exact cube supports, starting with all unpinned coordinates.

    The last entry is the residual; a repeated stationary entry is omitted.
    After t synchronous updates the deviation support is contained in trace[t].
    """
    pins = model.indices(pins)
    zero = model.relabeled()
    current = frozenset(range(model.n)) - pins
    trace = [current]
    while True:
        nxt = support_map(zero, current, pins)
        if not nxt <= current:
            raise AssertionError("Support descent was violated")
        if nxt == current:
            return tuple(trace)
        trace.append(nxt)
        current = nxt


def physical_residual_model(model, residual):
    """Exact dynamics on the invariant residual cube, reindexed 0,...,|R|-1.

    Returns (zero-target residual model or None for an empty residual,
    sorted original node indices). Coordinate relabeling precedes restriction.
    """
    residual = model.indices(residual)
    nodes = tuple(sorted(residual))
    if not nodes:
        return None, nodes
    locations = {node: index for index, node in enumerate(nodes)}
    zero = model.relabeled()
    rules = []
    for node in nodes:
        active, table = cube_row_table(zero, node, residual)
        rules.append(LocalRule(tuple(locations[p] for p in active), table))
    return BooleanModel(tuple(rules), (0,) * len(nodes)), nodes


def copy_parent(zero_model, row, support):
    """Identify a row that is exactly one coordinate on the restricted cube."""
    active, table = cube_row_table(zero_model, row, support)
    for j, parent in enumerate(active):
        if all(value == ((code >> j) & 1) for code, value in enumerate(table)):
            return parent
    return None


def functional_cycles(parent_map):
    """Canonical cycles of a total one-parent map (orientation is immaterial)."""
    if any(p not in parent_map for p in parent_map.values()):
        raise ValueError("A parent is outside the functional graph")
    finished, cycles = set(), []
    for start in sorted(parent_map):
        order, position = [], {}
        node = start
        while node not in finished and node not in position:
            position[node] = len(order)
            order.append(node)
            node = parent_map[node]
        if node in position:
            cycle = order[position[node]:]
            first = cycle.index(min(cycle))
            cycles.append(tuple(cycle[first:] + cycle[:first]))
        finished.update(order)
    return tuple(sorted(cycles))


def pure_copy_certificate(model, baseline_pins=()):
    """Exact conditional optimum when the unpinned residual consists of copies.

    This statement counts additional pins beyond the fixed baseline. It is an
    unconditional total-cost optimum only if the baseline is independently
    proved mandatory, as it is for identity rows in the supplied worked example.
    """
    baseline = model.indices(baseline_pins)
    trace = support_trace(model, baseline)
    residual = trace[-1]
    zero = model.relabeled()
    parents = {i: copy_parent(zero, i, residual) for i in sorted(residual)}
    if any(p is None for p in parents.values()):
        return None
    cycles = functional_cycles(parents)
    return {
        "baseline_pins": sorted(baseline),
        "support_trace": [sorted(s) for s in trace],
        "residual": sorted(residual),
        "copy_parent": {str(i): p for i, p in parents.items()},
        "cycles": [list(c) for c in cycles],
        "conditional_optimum": len(cycles),
        "one_optimal_additional_set": [c[0] for c in cycles],
    }


def identity_rows(model):
    """Rows equal to their own coordinate on every physical local assignment."""
    identities = []
    for i, rule in enumerate(model.rules):
        if i in rule.parents:
            j = rule.parents.index(i)
            if all(v == ((c >> j) & 1) for c, v in enumerate(rule.table)):
                identities.append(i)
    return tuple(identities)


def all_minimum_copy_sets(certificate):
    """All additional minimum-cardinality sets for a pure-copy residual."""
    return tuple(tuple(sorted(choice)) for choice in product(*certificate["cycles"]))


def load_json_model(document):
    """Load the explicit neutral JSON format shipped with this archive."""
    rules = tuple(LocalRule(tuple(row["parents_zero_based"]),
                            tuple(row["physical_truth_table"]))
                  for row in document["rules"])
    return BooleanModel(rules, tuple(document["target"]))
