"""Run the fixed Baseline protocol, keeping every draw and explicit cost bounds."""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import numpy as np
from engine import backend_info, certify, exact_residual_cost, generate_network, make_rng

FIELDS = ["mixture","n","replicate","residual_size","double_rows","certified",
          "exact_fallback","cost_lower","cost_upper","solver_status","solver_work","pins"]

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def one_cell(job):
    root, mixture, n, protocol = job
    root = Path(root)
    out = root / "results" / "raw" / f"{mixture}_n{n}.csv.gz"
    start = time.perf_counter()
    summary = {"mixture": mixture, "n": n, "draws": 0,
               "certificate_failures": 0, "exact_fallback": 0, "unresolved": 0}
    temporary = out.with_suffix(".partial")
    with temporary.open("wb") as binary:
        # Fix gzip timestamp and embedded filename for deterministic raw bytes.
        with gzip.GzipFile(filename="", mode="wb", fileobj=binary, mtime=0) as zipped:
            import io
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=FIELDS)
                writer.writeheader()
                for replicate in range(protocol["replicates_per_cell"]):
                    rng = make_rng(protocol["master_seed"],
                                   protocol["mixture_ids"][mixture], n, replicate)
                    p, b, d, q = generate_network(n, mixture, rng)
                    result = certify(p, b, d, backend="cpp")
                    fallback = 0
                    work = 0
                    if result["certified"]:
                        lo = hi = result["cost"]
                        pins = result["pins"]
                        status = "certified"
                    else:
                        summary["certificate_failures"] += 1
                        f = exact_residual_cost(
                            p, b, result["R"],
                            max_residual=protocol["uncertified_cases"]["maximum_residual_size"],
                            work_cap=protocol["uncertified_cases"]["maximum_work_units_per_draw"])
                        work = f["work"]
                        if f["exact"]:
                            lo = hi = f["cost"]
                            pins = f["pins"]
                            status = "exact"
                            fallback = 1
                            summary["exact_fallback"] += 1
                        else:
                            lo, hi = 0, result["residual_size"]
                            # Pinning every residual vertex is a valid feasible
                            # set giving the recorded upper bound.
                            pins = np.flatnonzero(result["R"])
                            status = f["status"]
                            summary["unresolved"] += 1
                    assert 0 <= lo <= hi <= result["residual_size"] <= n
                    assert len(pins) == hi
                    writer.writerow({
                        "mixture": mixture, "n": n, "replicate": replicate,
                        "residual_size": result["residual_size"],
                        "double_rows": result["double_rows"],
                        "certified": int(result["certified"]), "exact_fallback": fallback,
                        "cost_lower": lo, "cost_upper": hi, "solver_status": status,
                        "solver_work": work, "pins": ";".join(map(str,pins))})
                    summary["draws"] += 1
    os.replace(temporary, out)
    summary.update({"seconds": time.perf_counter()-start,
                    "sha256": sha(out), "file": str(out.relative_to(root))})
    return summary

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    root=args.root.resolve()
    protocol_path=root/"Protocol_Baseline.json"
    protocol=json.loads(protocol_path.read_text())
    validation=root/"results/engine_validation.json"
    if not validation.exists():
        raise RuntimeError("Independent validation must finish before the main run.")
    record=json.loads(validation.read_text())
    if record.get("status") != "passed":
        raise RuntimeError("Independent engine validation has not passed.")
    if backend_info()["backend"] != "cpp":
        raise RuntimeError("The fixed main workload requires the validated C++ backend.")
    raw=root/"results/raw"
    raw.mkdir(parents=True,exist_ok=True)
    if any(raw.glob("*.csv.gz")):
        raise RuntimeError("Refusing to overwrite an existing main run.")
    sources={p.name:sha(p) for p in (root/"simulation").iterdir()
             if p.suffix in (".py",".cpp")}
    environment={"python":sys.version, "numpy":np.__version__, "platform":platform.platform(),
                 "compiler":subprocess.check_output(["g++","--version"],text=True).splitlines()[0],
                 "cpu_count":os.cpu_count(), "workers":args.workers}
    metadata={"protocol_sha256":sha(protocol_path), "sources_sha256":sources,
              "engine_validation_sha256":sha(validation),
              "environment":environment,"backend":backend_info(),
              "start_time_unix":time.time(),"status":"running","cells":[]}
    manifest=root/"results/run_manifest.json"
    manifest.write_text(json.dumps(metadata,indent=2)+"\n")
    jobs=[(str(root),m,n,protocol) for m in protocol["mixtures"] for n in protocol["sizes"]]
    start=time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(one_cell,j):j[1:3] for j in jobs}
        for future in as_completed(futures):
            cell=future.result()
            metadata["cells"].append(cell)
            manifest.write_text(json.dumps(metadata,indent=2)+"\n")
            print(json.dumps({"completed":len(metadata["cells"]),"of":len(jobs),**cell}),flush=True)
    assert sum(c["draws"] for c in metadata["cells"])==protocol["total_networks"]
    assert sha(protocol_path)==metadata["protocol_sha256"]
    # Analysis/export scripts may be developed while the main run is active;
    # random generation, physical evaluation and the runner remain frozen.
    for name in ["engine.py","engine_core.cpp","run_experiments.py"]:
        assert sha(root/"simulation"/name)==sources[name],f"Main-run source changed: {name}"
    metadata.update({"status":"completed","wall_seconds":time.perf_counter()-start,
                     "end_time_unix":time.time()})
    metadata["cells"].sort(key=lambda c:(c["mixture"],c["n"]))
    manifest.write_text(json.dumps(metadata,indent=2)+"\n")
    print(json.dumps({"status":"completed","draws":protocol["total_networks"],
                      "seconds":metadata["wall_seconds"]}),flush=True)

if __name__=="__main__":
    main()
