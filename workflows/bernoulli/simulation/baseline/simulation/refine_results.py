"""Replay only previously unresolved draws and tighten deterministic intervals.

The original raw data and fixed protocol are immutable. This analysis-only
refinement is separately declared in Bound_Refinement_Protocol.json.
"""
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import time
from engine import make_rng, generate_network
from refine_control import control_bounds

ROOT=Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=json.loads((ROOT/"Protocol_Baseline.json").read_text())
    validation=json.loads((ROOT/"results/refined_bounds_validation.json").read_text())
    assert validation["status"]=="passed"
    directory=ROOT/"results/refined/raw"
    directory.mkdir(parents=True,exist_ok=True)
    if any(directory.glob("*.csv.gz")):
        raise RuntimeError("Refusing to overwrite an existing refined analysis.")
    manifest={"status":"running","new_random_draws":0,"replayed_unresolved_draws":0,
              "protocol_sha256":sha(ROOT/"Protocol_Baseline.json"),
              "analysis_addendum_sha256":sha(ROOT/"Bound_Refinement_Protocol.json"),
              "source_sha256":{n:sha(ROOT/"simulation"/n) for n in
                               ["refine_control.py","refine_results.py","engine.py","engine_core.cpp"]},
              "validation_sha256":sha(ROOT/"results/refined_bounds_validation.json"),
              "files":[]}
    start=time.perf_counter()
    total=0
    for source in sorted((ROOT/"results/raw").glob("*.csv.gz")):
        out=directory/source.name
        count=0
        before_hash=sha(source)
        with gzip.open(source,"rt",newline="") as f, out.open("wb") as out_binary:
            reader=csv.DictReader(f)
            with gzip.GzipFile(filename="",mode="wb",fileobj=out_binary,mtime=0) as gz:
                with io.TextIOWrapper(gz,encoding="utf-8",newline="") as dest:
                    writer=csv.DictWriter(dest,fieldnames=reader.fieldnames)
                    writer.writeheader()
                    for row in reader:
                        if row["cost_lower"] != row["cost_upper"]:
                            n=int(row["n"]);m=row["mixture"];rep=int(row["replicate"])
                            rng=make_rng(p["master_seed"],p["mixture_ids"][m],n,rep)
                            parents,rules,degrees,q=generate_network(n,m,rng)
                            b=control_bounds(parents,rules,degrees)
                            assert b["residual_size"]==int(row["residual_size"])
                            assert b["double_rows"]==int(row["double_rows"])
                            assert int(row["cost_lower"])<=b["lower"]<=b["upper"]<=int(row["cost_upper"])
                            assert len(b["pins"])==b["upper"]
                            row.update(cost_lower=b["lower"],cost_upper=b["upper"],
                                       pins=";".join(map(str,b["pins"])),
                                       solver_status="refined_"+row["solver_status"])
                            count+=1
                        writer.writerow(row)
                        total+=1
        assert sha(source)==before_hash
        manifest["replayed_unresolved_draws"]+=count
        manifest["files"].append({"file":source.name,"original_sha256":before_hash,
                                  "refined_sha256":sha(out),"refined_draws":count})
        print(json.dumps({"file":source.name,"replayed":count}),flush=True)
    assert total==p["total_networks"] and manifest["replayed_unresolved_draws"]==16920
    manifest.update(status="completed",total_draws=total,seconds=time.perf_counter()-start)
    (ROOT/"results/refined/refinement_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps({"status":"completed","total_draws":total,
                      "replayed_unresolved_draws":manifest["replayed_unresolved_draws"]}),flush=True)

if __name__=="__main__":
    main()
