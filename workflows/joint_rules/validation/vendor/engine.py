"""Reproducible planted Boolean configuration networks and exact certification.

The random generator is NumPy PCG64 supplied by the caller.  The C++ accelerator
only evaluates deterministic support/cycle calculations.  Set backend='python'
to run the equivalent linear-work reference implementation without a compiler.
The exact fallback enumerates physical states, not the conservative slot cube.
"""
from __future__ import annotations

import ctypes
import hashlib
import itertools
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np

MIXTURES = {
    "A": ((1, 0.50, 0.5), (3, 0.10, 0.5)),
    "B": ((1, 0.10, 0.5), (3, 0.50, 0.5)),
    "C": ((1, 0.80, 0.5), (3, 0.20, 0.5)),
    "D": ((1, 0.35, 0.5), (3, 0.35, 0.5)),
    "E": ((0, 0.80, 0.5), (4, 0.30, 0.5)),
}
MIXTURE_IDS = {name: i for i, name in enumerate(MIXTURES, start=1)}
_CORE = None
_CORE_ERROR = None


def make_rng(master_seed: int, mixture_id: int, n: int, replicate: int):
    """Each draw is independent of worker count, task order, and prior draws."""
    return np.random.Generator(np.random.PCG64(
        np.random.SeedSequence([master_seed, mixture_id, n, replicate])))


def generate_network(n, mixture, rng):
    """Generate a uniform matching with fixed vertex types and independent rules.

    ``mixture`` is A--E or a sequence of (output_degree, rule_probability,
    proportion).  Type counts must be exact integers and output stubs total 2n.
    The type order determines fixed vertex labels.  A uniform permutation of
    the repeated-source array is the projection of a uniform permutation of
    distinct labeled output stubs.  No simplicity rejection is performed.

    RNG call order: rng.permutation(output_sources), then rng.random((n, 3)).
    Rule columns are g(1,0), g(0,1), g(1,1); g(0,0) is identically zero.
    """
    if isinstance(n, bool) or int(n) != n or n < 1:
        raise ValueError("n must be a positive integer")
    n = int(n)
    types = MIXTURES[mixture] if isinstance(mixture, str) else tuple(mixture)
    if not types:
        raise ValueError("at least one type is required")
    counts, ds, qs = [], [], []
    for d, q, fraction in types:
        if int(d) != d or d < 0 or not (0 < q < 1) or fraction <= 0:
            raise ValueError("invalid degree, probability, or proportion")
        k = int(round(n * fraction))
        if not np.isclose(k, n * fraction, atol=1e-10, rtol=0):
            raise ValueError("n does not permit the specified exact type counts")
        counts.append(k)
        ds.append(int(d))
        qs.append(float(q))
    if sum(counts) != n or sum(k*d for k, d in zip(counts, ds)) != 2*n:
        raise ValueError("type counts must total n and output stubs must total 2n")
    degrees = np.repeat(np.asarray(ds, dtype=np.int64), counts)
    q = np.repeat(np.asarray(qs, dtype=np.float64), counts)
    output_sources = np.repeat(np.arange(n, dtype=np.int64), degrees)
    parents = np.asarray(rng.permutation(output_sources), dtype=np.int64).reshape(n, 2)
    rules = np.asarray(rng.random((n, 3)) < q[:, None], dtype=np.uint8)
    return parents, rules, degrees, q


def _arrays(parents, rules, degrees=None):
    parents_raw, rules_raw = np.asarray(parents), np.asarray(rules)
    if parents_raw.ndim != 2 or parents_raw.shape[1] != 2:
        raise ValueError("parents must have shape (n,2)")
    n = parents_raw.shape[0]
    if rules_raw.shape != (n, 3):
        raise ValueError("rules must have shape (n,3)")
    if np.any((rules_raw != 0) & (rules_raw != 1)):
        raise ValueError("rules must be binary")
    if np.any(parents_raw != parents_raw.astype(np.int64)):
        raise ValueError("parent indices must be integers")
    p = np.ascontiguousarray(parents_raw, dtype=np.int64)
    b = np.ascontiguousarray(rules_raw, dtype=np.uint8)
    if np.any(p < 0) or np.any(p >= n):
        raise ValueError("parent index outside [0,n)")
    if degrees is not None:
        d = np.asarray(degrees)
        if d.shape != (n,) or not np.array_equal(np.bincount(p.ravel(), minlength=n), d):
            raise ValueError("degrees do not equal the realized output multiplicities")
    return p, b


def build_core():
    """Compile the bundled deterministic accelerator; output is source-hashed."""
    source = Path(__file__).with_name("engine_core.cpp")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()[:16]
    build_dir = source.parent / "_build"
    build_dir.mkdir(exist_ok=True)
    library = build_dir / ("engine_core_" + digest + ".so")
    if not library.exists():
        fd, tmp = tempfile.mkstemp(prefix="engine_core_", suffix=".so", dir=build_dir)
        os.close(fd)
        try:
            result = subprocess.run(
                ["g++", "-std=c++17", "-O3", "-fPIC", "-shared", str(source), "-o", tmp],
                capture_output=True, text=True, check=False)
            if result.returncode:
                raise RuntimeError("C++ build failed: " + result.stderr)
            os.replace(tmp, library)
        finally:
            if os.path.exists(tmp):
                os.unlink(tmp)
    return library


def _load_core():
    global _CORE, _CORE_ERROR
    if _CORE is None and _CORE_ERROR is None:
        try:
            lib = ctypes.CDLL(str(build_core()))
            fn = lib.certify_core
            p64 = np.ctypeslib.ndpointer(dtype=np.int64, flags="C_CONTIGUOUS")
            p8 = np.ctypeslib.ndpointer(dtype=np.uint8, flags="C_CONTIGUOUS")
            fn.argtypes = [ctypes.c_int64, p64, p8, p8, p64, p64, p64]
            fn.restype = ctypes.c_int
            _CORE = fn
        except (OSError, RuntimeError) as exc:
            _CORE_ERROR = str(exc)
    return _CORE


def backend_info():
    core = _load_core()
    return {"backend": "cpp" if core is not None else "python",
            "compiler_error": _CORE_ERROR,
            "cpp_sha256": hashlib.sha256(Path(__file__).with_name("engine_core.cpp").read_bytes()).hexdigest()}


def _certify_python(p, b):
    n = len(p)
    flat = p.ravel().tolist()
    rules = b.tolist()
    counts = [0] * n
    for source in flat:
        counts[source] += 1
    offset = [0] * (n+1)
    for i in range(n):
        offset[i+1] = offset[i] + counts[i]
    cursor = offset[:-1].copy()
    slots = [0] * (2*n)
    for e, source in enumerate(flat):
        slots[cursor[source]] = e
        cursor[source] += 1
    alive = [bool(a or c or d) for a, c, d in rules]
    mask = [3] * n
    queue = [i for i in range(n) if not alive[i]]
    head = 0
    while head < len(queue):
        source = queue[head]
        head += 1
        for j in range(offset[source], offset[source+1]):
            e = slots[j]
            i, s = divmod(e, 2)
            if not alive[i]:
                continue
            mask[i] &= ~(1 << s)
            m = mask[i]
            if not (m == 3 or (m == 1 and rules[i][0]) or (m == 2 and rules[i][1])):
                alive[i] = False
                queue.append(i)
    residual_parent = [-1] * n
    doubles = 0
    for i in range(n):
        if not alive[i]:
            continue
        first, second = alive[flat[2*i]], alive[flat[2*i+1]]
        if first and second:
            doubles += 1
        elif first:
            residual_parent[i] = flat[2*i]
        elif second:
            residual_parent[i] = flat[2*i+1]
        else:
            raise AssertionError("support fixed point has a row without residual slots")
    chosen = np.zeros(n, dtype=bool)
    if not doubles:
        seen = [0] * n
        for i in range(n):
            if not alive[i] or seen[i]:
                continue
            stamp, v = i+1, i
            while not seen[v]:
                seen[v] = stamp
                v = residual_parent[v]
            if seen[v] == stamp:
                entry, smallest = v, v
                v = residual_parent[v]
                while v != entry:
                    smallest = min(smallest, v)
                    v = residual_parent[v]
                chosen[smallest] = True
    pins = np.flatnonzero(chosen)
    return {"R": np.asarray(alive, dtype=bool), "residual_size": sum(alive),
            "double_rows": doubles, "certified": doubles == 0,
            "pins": pins,
            "cost": len(pins) if not doubles else -1,
            "residual_parent": np.asarray(residual_parent, dtype=np.int64)}


def certify(parents, rules, degrees=None, backend="auto"):
    """Greatest conservative slot support and the exact copy-cycle certificate.

    Both implementations perform O(n) deterministic work and count residual
    input slots, including both slots when they share a physical source.
    A failed certificate is an unknown optimum, not a failed intervention.
    """
    if backend not in ("auto", "cpp", "python"):
        raise ValueError("backend must be auto, cpp, or python")
    p, b = _arrays(parents, rules, degrees)
    core = _load_core() if backend != "python" else None
    if core is None:
        if backend == "cpp":
            raise RuntimeError(_CORE_ERROR or "C++ backend unavailable")
        return _certify_python(p, b)
    n = len(p)
    alive = np.empty(n, dtype=np.uint8)
    residual_parent = np.empty(n, dtype=np.int64)
    pins = np.empty(n, dtype=np.int64)
    stats = np.empty(4, dtype=np.int64)
    code = core(n, p, b, alive, residual_parent, pins, stats)
    if code:
        raise RuntimeError("C++ certificate failed with internal status " + str(code))
    r, doubles, cost, npins = map(int, stats)
    return {"R": alive.astype(bool), "residual_size": r, "double_rows": doubles,
            "certified": doubles == 0, "pins": pins[:npins].copy(), "cost": cost,
            "residual_parent": residual_parent}


def exact_residual_cost(parents, rules, R, max_residual=12, work_cap=2_000_000):
    """Exact minimum on the residual cube, with deterministic resource bounds.

    R must be the greatest slot support returned by certify.  Outside-R vertices
    remain zero on this cube.  Uniform capture makes the residual optimum equal
    to the full physical optimum, including when the certificate fails.
    The returned status distinguishes resource limits from an exact result.

    Work units are residual-row evaluations when constructing every physical
    successor, plus controlled successor-edge traversals during global tests.
    Compatible states have every pinned coordinate zero at time zero.
    """
    p, b = _arrays(parents, rules)
    R = np.asarray(R, dtype=bool)
    n = len(p)
    if R.shape != (n,):
        raise ValueError("R must have shape (n,)")
    if int(max_residual) != max_residual or max_residual < 0:
        raise ValueError("max_residual must be a nonnegative integer")
    if int(work_cap) != work_cap or work_cap < 0:
        raise ValueError("work_cap must be a nonnegative integer")
    work_cap = int(work_cap)
    vertices = np.flatnonzero(R)
    r = len(vertices)
    empty = np.empty(0, dtype=np.int64)
    def unresolved(status, work):
        return {"exact": False, "cost": -1, "pins": empty.copy(),
                "status": status, "work": int(work)}
    if r == 0:
        return {"exact": True, "cost": 0, "pins": empty,
                "status": "exact", "work": 0}
    if r > max_residual:
        return unresolved("residual_cap", 0)
    count = 1 << r
    initial_work = r * count
    if initial_work > work_cap or r >= 63:
        return unresolved("work_cap", 0)
    local = np.full(n, -1, dtype=np.int64)
    local[vertices] = np.arange(r)
    states = np.arange(count, dtype=np.int64)
    successor = np.zeros(count, dtype=np.int64)
    for j, i in enumerate(vertices):
        u, v = local[p[i, 0]], local[p[i, 1]]
        x = (states >> u) & 1 if u >= 0 else np.zeros(count, dtype=np.int64)
        y = (states >> v) & 1 if v >= 0 else np.zeros(count, dtype=np.int64)
        # Repeated parents use the same physical coordinate in both expressions.
        output = (((x == 1) & (y == 0) & bool(b[i, 0])) |
                  ((x == 0) & (y == 1) & bool(b[i, 1])) |
                  ((x == 1) & (y == 1) & bool(b[i, 2])))
        successor |= output.astype(np.int64) << j
    work = initial_work
    successor = successor.tolist()
    for k in range(r+1):
        for selected in itertools.combinations(range(r), k):
            pinmask = sum(1 << j for j in selected)
            free_mask = (count-1) ^ pinmask
            color = bytearray(count)
            color[0] = 2
            success = True
            for start in range(1, count):
                if start & pinmask or color[start]:
                    continue
                path = []
                x = start
                while color[x] == 0:
                    if work == work_cap:
                        return unresolved("work_cap", work)
                    color[x] = 1
                    path.append(x)
                    x = successor[x] & free_mask
                    work += 1
                if color[x] == 1:
                    success = False
                    break
                for x in path:
                    color[x] = 2
            if success:
                return {"exact": True, "cost": k,
                        "pins": vertices[np.asarray(selected, dtype=np.int64)].copy(),
                        "status": "exact", "work": int(work)}
    raise AssertionError("pinning all residual vertices must succeed")


if __name__ == "__main__":
    print(backend_info())
