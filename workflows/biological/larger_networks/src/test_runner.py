"""Small-state independent checks of search logic and declared cap behavior."""
import itertools
import json
import time

from benchmark import (BooleanModel, LocalRule, WorkBudget, classify_pins,
                       discover_target, exact_search, parse_boolean)


def direct_orbit_oracle(successor, pins):
    mask = sum(1 << p for p in pins)
    maximum, success = 0, True
    for initial in range(len(successor)):
        if initial & mask:
            continue
        state, visited, steps = initial, set(), 0
        while state != 0 and state not in visited:
            visited.add(state)
            state = successor[state] & ~mask
            steps += 1
        if state:
            success = False
        else:
            maximum = max(maximum, steps)
    return success, maximum


def run_tests():
    checked = 0
    for tail in itertools.product(range(4), repeat=3):
        successor = [0, *tail]
        for cardinality in range(3):
            for pins in itertools.combinations(range(2), cardinality):
                expected, hitting = direct_orbit_oracle(successor, pins)
                obtained = classify_pins(successor, 2, pins,
                                         WorkBudget(100, time.perf_counter() + 5))
                assert obtained['globally_convergent'] == expected
                if expected:
                    assert obtained['max_residual_hitting_time'] == hitting
                else:
                    cycle = obtained['failure_cycle']
                    mask = sum(1 << p for p in pins)
                    assert cycle and 0 not in cycle
                    assert all(successor[a] & ~mask == b for a, b in zip(cycle, cycle[1:] + cycle[:1]))
                checked += 1
    budgets = {'max_residual_nodes': 18, 'max_transition_evaluations': 10000,
               'max_additional_cost': 6}
    # Successful empty pin set with a nonempty residual; additional pins can fail.
    nonlinear = BooleanModel((LocalRule((0, 1), (0, 1, 0, 0)),
                              LocalRule((0,), (0, 1))), (0, 0))
    report = exact_search(nonlinear, (), ('a', 'b'), budgets, time.perf_counter() + 5)
    assert report['residual_size'] == 2
    assert report['search']['lower_bound'] == report['search']['upper_bound'] == 0
    assert report['search']['all_minimum_additional_sets'] == [[]]
    # No budget can silently convert a feasible upper bound into an exact optimum.
    two_identities = BooleanModel((LocalRule((0,), (0, 1)), LocalRule((1,), (0, 1))), (0, 0))
    cap_cases = []
    for changes, deadline, expected, bounds in (
        ({'max_transition_evaluations': 1}, time.perf_counter() + 5, 'transition_budget', [0, 2]),
        ({'max_additional_cost': 0}, time.perf_counter() + 5, 'additional_cost_cap', [1, 2]),
        ({}, time.perf_counter() - 1, 'wall_budget', [0, 2]),
    ):
        report = exact_search(two_identities, (), ('a', 'b'), {**budgets, **changes}, deadline)
        search = report['search']
        assert search['status'] == expected
        assert [search['lower_bound'], search['upper_bound']] == bounds
        assert not search['all_minimum_sets_complete']
        cap_cases.append(expected)
    # Pure-copy large-residual fallback supplies a mathematical exact cost.
    report = exact_search(two_identities, (), ('a', 'b'), {**budgets, 'max_residual_nodes': 1}, time.perf_counter() + 5)
    assert report['search']['status'] == 'exact_copy_certificate'
    assert report['search']['lower_bound'] == report['search']['upper_bound'] == 2
    # Fixed-target discovery may fail at a periodic orbit without proving no fixed point.
    discovery = discover_target((LocalRule((0,), (1, 0)),), ('a',), (), {}, 8, time.perf_counter() + 5)
    assert discovery['status'] == 'target_not_found_by_declared_orbit'
    assert discovery['cycle_length'] == 2 and discovery['target'] is None
    rejected = 0
    for formula in ('__import__("os")', 'a[0]', 'a + b', '3'):
        try:
            parse_boolean(formula)
        except (ValueError, SyntaxError):
            rejected += 1
    assert rejected == 4
    return {'status': 'PASS', 'all_two_coordinate_zero_fixed_map_pin_cases': checked,
            'nonmonotone_empty_pin_success_preserved': True,
            'cap_cases_checked': cap_cases, 'copy_fallback_exact_cost_checked': True,
            'periodic_orbit_discovery_status_checked': True,
            'unsafe_or_nonboolean_parser_inputs_rejected': rejected}


if __name__ == '__main__':
    print(json.dumps(run_tests(), indent=2))
