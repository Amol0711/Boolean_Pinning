"""Deterministic control-cost bounds from the original residual double rows.

This is a separate post-run refinement.  It does not alter the random generator,
the original certificate, the capped solver, or the frozen main-run records.
"""
from __future__ import annotations

import numpy as np

from engine import certify


def control_bounds(parents, rules, degrees=None):
    """Return rigorous bounds and a successful pin set attaining the upper bound.

    Let R be the greatest original slot support and D its vertices with two
    residual slots.  Replace each row in D by the constant-zero rule and compute
    its residual copy-cycle count c.  Then c <= CK <= c + |D|.  The upper pin set
    contains D and one vertex from each surviving copy cycle.  For D empty this
    recovers the original exact certificate.  Source arrays are never mutated.

    The lower bound uses disjoint cycles consisting of original unary residual
    rows.  Each remains physically invariant at one on the original R cube
    unless that cycle itself is pinned.  The upper bound is certified by the
    constant-zero replacement followed by one pin per remaining copy cycle.
    """
    original = certify(parents, rules, degrees)
    p = np.asarray(parents, dtype=np.int64)
    R = original["R"]
    double_mask = R & R[p[:, 0]] & R[p[:, 1]]
    double_vertices = np.flatnonzero(double_mask)
    assert len(double_vertices) == original["double_rows"]
    if original["certified"]:
        assert len(double_vertices) == 0
        c = original["cost"]
        return {"lower": c, "upper": c, "exact": True,
                "pins": original["pins"].copy(),
                "original_certified": True,
                "residual_size": original["residual_size"],
                "double_rows": 0, "double_vertices": double_vertices,
                "cycle_cost": c,
                "refined_residual_size": original["residual_size"]}
    changed_rules = np.array(rules, dtype=np.uint8, copy=True)
    changed_rules[double_mask] = 0
    after = certify(p, changed_rules, degrees)
    # Monotonicity of support under row removal makes every survivor an original
    # unary residual row; constant-zero rows cannot survive in the new support.
    assert np.all(~after["R"] | R)
    assert not np.any(after["R"] & double_mask)
    assert after["certified"] and after["double_rows"] == 0
    assert not np.any(double_mask[after["pins"]])
    c = after["cost"]
    pin_mask = double_mask.copy()
    pin_mask[after["pins"]] = True
    pins = np.flatnonzero(pin_mask)
    lower, upper = c, c + len(double_vertices)
    assert len(pins) == upper
    assert 0 <= lower <= upper <= original["residual_size"]
    return {"lower": lower, "upper": upper, "exact": lower == upper,
            "pins": pins, "original_certified": False,
            "residual_size": original["residual_size"],
            "double_rows": len(double_vertices),
            "double_vertices": double_vertices, "cycle_cost": c,
            "refined_residual_size": after["residual_size"]}
