"""Build the two prespecified larger-network panels from frozen certificates.

Panel A is the exact coordinate-cube support envelope for CD4_APC.
Panel B displays every varying coordinate of the MAC_BACTERIA four-cycle.
Its full trajectory and all four cycle edges are independently replayed from
native BNET rules. The two panels have different, explicitly labeled semantics.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
PANELS = (
    ("CD4_APC", "A. CD4 (188 nodes)", "APC and MHC II active"),
    ("MAC_BACTERIA", "B. Macrophage (321 nodes)", "Bacteria and External Activator active"),
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tex_escape(text):
    return str(text).replace("_", r"\_").replace("&", r"\&")


def load_panel(results, context_id, title, context_label, data_dir):
    path = results / context_id / "certificate.json"
    certificate = json.loads(path.read_text())
    if certificate["context_id"] != context_id:
        raise ValueError("Certificate context ID differs from its path")
    metadata = {
        "context_id": context_id,
        "source_certificate": str(path.relative_to(ROOT)),
        "source_certificate_sha256": digest(path),
        "node_count": len(certificate["node_order"]),
        "input_count": len(certificate["input_node_names"]),
        "input_assignment": certificate["baseline_input_assignment"],
        "target_discovery_status": certificate["discovery"]["status"],
        "ordinate_semantics": "Exact coordinate-cube support cardinality; an outer bound on possible deviations, not an individual state trajectory or exact reachable-support statistic.",
    }
    if certificate["discovery"]["target"] is None:
        return load_cycle_panel(certificate, metadata, title, context_label, data_dir)

    baseline = certificate["baseline_support_sizes"]
    if baseline != [len(s) for s in certificate["baseline_support_trace"]]:
        raise ValueError("Baseline support cardinalities do not match their sets")
    expected_start = metadata["node_count"] - metadata["input_count"]
    if not baseline or baseline[0] != expected_start:
        raise ValueError("Initial support must exclude all clamped inputs")
    search = certificate["search"]
    exact = search["lower_bound"] == search["upper_bound"]
    optimal_trace = certificate.get("first_feasible_full_support_trace")
    optimum = [len(s) for s in optimal_trace] if exact and optimal_trace is not None else None
    if optimum is not None:
        if optimum[0] != expected_start - search["upper_bound"]:
            raise ValueError("Initial optimal support does not match the exact pin count")
    for values in (baseline, optimum):
        if values is not None and any(b > a for a, b in zip(values, values[1:])):
            raise ValueError("Coordinate-cube support must descend")

    same = optimum == baseline
    last_round = max(len(baseline), len(optimum or [])) - 1
    # Include one unchanged round to make a nonempty stationary limit explicit.
    display_last = last_round + 1
    csv_path = data_dir / ("fig5_" + context_id.lower() + ".csv")
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("support_round", "baseline_support", "optimal_support",
                         "baseline_stationary_extension", "optimal_stationary_extension"))
        for k in range(display_last + 1):
            writer.writerow((k, baseline[min(k, len(baseline) - 1)],
                             optimum[min(k, len(optimum) - 1)] if optimum is not None else "",
                             int(k >= len(baseline)),
                             int(k >= len(optimum)) if optimum is not None else ""))
    metadata.update({
        "plotted": True,
        "csv_file": str(csv_path.relative_to(ROOT)),
        "csv_sha256": digest(csv_path),
        "exact_additional_cost": search["upper_bound"] if exact else None,
        "additional_cost_bounds": [search["lower_bound"], search["upper_bound"]],
        "search_status": search["status"],
        "selected_additional_pin_values": search["feasible_additional_values"],
        "baseline_original_sizes": baseline,
        "optimal_original_sizes": optimum,
        "baseline_equals_optimal": same,
        "baseline_residual_size": baseline[-1],
        "optimal_support_limit": optimum[-1] if optimum is not None else None,
        "display_last_round": display_last,
        "stationary_extension_note": "A support is repeated only after its exact support-map trace has reached a fixed point.",
    })
    return {
        "title": title, "context_label": context_label, "metadata": metadata,
        "csv": "data/" + csv_path.name, "baseline": baseline,
        "optimum": optimum, "same": same, "last": display_last,
        "cost": search["upper_bound"] if exact else None,
        "ymax": int(math.ceil(baseline[0] / 20) * 20),
    }


def load_cycle_panel(certificate, metadata, title, context_label, data_dir):
    """Check source dynamics independently and export the actual cycle."""
    sys.path.insert(0, str(ROOT / "larger_networks/verification"))
    from independent_core import read_bnet, step

    if certificate["context_id"] != "MAC_BACTERIA":
        raise ValueError("Cycle panel is reserved for the frozen macrophage context")
    source = ROOT / "larger_networks/data/macrophage_321/model.bnet"
    if digest(source) != certificate["source"]["model_bnet_sha256"]:
        raise ValueError("Native BNET source hash differs from the certificate")
    rows, regulated, inputs = read_bnet(source)
    order = certificate["node_order"]
    if order != regulated + inputs or inputs != certificate["input_node_names"]:
        raise ValueError("Native BNET coordinate order differs from the certificate")
    discovery = certificate["discovery"]
    initial = {name: 0 for name in regulated}
    initial.update(certificate["baseline_input_assignment"])
    if [initial[name] for name in order] != discovery["initial_state"]:
        raise ValueError("Declared initial state differs from the frozen protocol")
    # Replay until the first repeated full state, without importing the engine.
    seen, trajectory = {}, []
    state = initial
    while tuple(state[name] for name in order) not in seen:
        bits = tuple(state[name] for name in order)
        seen[bits] = len(trajectory)
        trajectory.append(list(bits))
        if len(trajectory) > certificate["budget"]["target_orbit_max_updates"]:
            raise ValueError("Independent orbit exceeded the declared target budget")
        state = step(rows, state)
    repeated = [state[name] for name in order]
    cycle_start = seen[tuple(repeated)]
    cycle = trajectory[cycle_start:]
    if (cycle != discovery["cycle"] or cycle_start != discovery["cycle_start_update"]
            or len(cycle) != discovery["cycle_length"]
            or len(trajectory) != discovery["updates"]):
        raise ValueError("Independent full-state orbit disagrees with stored cycle")
    if len(cycle) != 4 or cycle_start != 11 or len(set(map(tuple, cycle))) != 4:
        raise ValueError("Frozen four-cycle data changed")
    trajectory.append(repeated)
    varying = [i for i in range(len(order)) if len({state[i] for state in cycle}) > 1]
    constants = [i for i in range(len(order)) if i not in varying]
    if len(varying) != 19 or len(constants) != 302:
        raise ValueError("Unexpected varying/constant coordinate counts")
    if any(order[i] in inputs for i in varying):
        raise ValueError("An environmental input varies on the recorded cycle")
    phases = list(range(cycle_start, cycle_start + len(cycle) + 1))
    displayed = cycle + [cycle[0]]
    for phase, value in zip(phases, displayed):
        if trajectory[phase] != value:
            raise ValueError("Displayed cycle phase differs from independent orbit")
    orbit_path = data_dir / "fig5_mac_bacteria_full_orbit.csv"
    with orbit_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["synchronous_update"] + order)
        writer.writerows([[t] + value for t, value in enumerate(trajectory)])
    matrix_path = data_dir / "fig5_mac_bacteria_cycle.csv"
    with matrix_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["node_index_zero_based", "source_node_name", "display_label"] +
                        ["update_" + str(t) for t in phases])
        writer.writerows([[i, order[i], order[i].removeprefix("v_")] +
                          [state[i] for state in displayed] for i in varying])
    constants_path = data_dir / "fig5_mac_bacteria_constant_coordinates.csv"
    with constants_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["node_index_zero_based", "source_node_name", "value_on_cycle",
                         "external_input"])
        writer.writerows([[i, order[i], cycle[0][i], int(order[i] in inputs)] for i in constants])
    # Keep the status CSV for backwards-compatible provenance, never as plot data.
    status_path = data_dir / "fig5_mac_bacteria_status.csv"
    with status_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("context_id", "target_discovery_status", "cycle_start_update",
                         "cycle_length", "fixed_target_certificate_applied"))
        writer.writerow((certificate["context_id"], discovery["status"], cycle_start,
                         len(cycle), False))
    metadata.update({
        "plotted": True, "panel_kind": "verified_binary_cycle_raster",
        "ordinate_semantics": "Named coordinates; binary states of one native-source orbit, not target-relative support or a control certificate.",
        "abscissa_semantics": "Actual synchronous update number from the declared zero-regulated initial state with frozen input values.",
        "native_source_sha256": digest(source),
        "independent_parser_sha256": digest(ROOT / "larger_networks/verification/independent_core.py"),
        "independent_replay": "PASS; all 15 transitions from update 0 through update 15 checked, including all four cycle edges",
        "declared_orbit_cycle_start_update": cycle_start,
        "declared_orbit_cycle_length": len(cycle),
        "displayed_updates": phases, "repeated_first_phase_update": phases[-1],
        "varying_coordinate_count": len(varying),
        "varying_coordinate_names": [order[i] for i in varying],
        "constant_coordinate_count": len(constants),
        "constant_regulated_coordinate_count": len(constants) - len(inputs),
        "constant_input_count": len(inputs),
        "constant_zero_count": sum(cycle[0][i] == 0 for i in constants),
        "constant_one_count": sum(cycle[0][i] == 1 for i in constants),
        "label_transformation": "Only the common source prefix v_ is omitted.",
        "fixed_target_certificate_applied": False,
        "csv_files": {str(path.relative_to(ROOT)): digest(path)
                      for path in (orbit_path, matrix_path, constants_path, status_path)},
        "cycle_values_by_name": {order[i]: [value[i] for value in displayed] for i in varying},
    })
    return {"title": title, "context_label": context_label, "metadata": metadata,
            "kind": "cycle", "phases": phases,
            "labels": [order[i].removeprefix("v_") for i in varying],
            "values": [[value[i] for value in displayed] for i in varying]}


def support_panel_tex(panel):
    last = panel["last"]
    result = r"""\begin{axis}[name=panelA,at={(0,0)},anchor=south west,
width=6.2cm,height=6.4cm,scale only axis,
title={A. CD4 (188 nodes)},
xmin=-0.2,xmax=XMAX,ymin=-4,ymax=YMAX,
xtick={0,2,...,XTICKMAX},ytick={0,40,...,160},
xlabel={Support-reduction round $k$},
ylabel={Potentially non-target coordinates},
grid=major,grid style={black!10},
legend style={font=\scriptsize,draw=none,fill=white,
at={(0.97,0.97)},anchor=north east},legend cell align=left]
\addplot[curveorange,thick,dashed,mark=square*,mark size=1.8pt]
table[x=support_round,y=baseline_support,col sep=comma]{CSV};
\addlegendentry{Inputs clamped}
\addplot[curveblue,thick,mark=*,mark size=1.8pt]
table[x=support_round,y=optimal_support,col sep=comma]{CSV};
\addlegendentry{Inputs + minimum pins}
\node[anchor=north east,align=right,font=\footnotesize,text=black!80]
at (axis description cs:0.96,0.81)
{$|B|=34$, $c_B=4$\\Baseline residual $|R_B|=17$};
\end{axis}
""".replace("XMAX", str(last + 0.25)).replace("YMAX", str(panel["ymax"] + 8)).replace(
        "XTICKMAX", str(last)).replace("CSV", panel["csv"])
    if panel["cost"] != 4 or panel["metadata"]["baseline_residual_size"] != 17:
        raise ValueError("Panel A annotation no longer matches its certificate")
    return result


def cycle_panel_tex(panel):
    nrows = len(panel["labels"])
    labels = ",".join("{" + tex_escape(label) + "}" for label in panel["labels"])
    ticks = ",".join(str(i) for i in range(nrows))
    result = r"""\begin{axis}[name=panelB,at={(10.6cm,0)},anchor=south west,
width=4.5cm,height=6.4cm,scale only axis,
title={B. Macrophage (321 nodes)},
xmin=10.5,xmax=15.5,ymin=-0.5,ymax=18.5,y dir=reverse,
xtick={11,12,13,14,15},ytick={TICKS},yticklabels={LABELS},
yticklabel style={font=\scriptsize},
xticklabel style={font=\footnotesize},
xlabel={Synchronous update $t$},
axis on top,tick align=outside,clip=false]
""".replace("TICKS", ticks).replace("LABELS", labels)
    for row, values in enumerate(panel["values"]):
        for phase, value in zip(panel["phases"], values):
            color = "curveblue" if value else "white"
            ink = "white" if value else "black!75"
            result += (r"\path[fill=" + color + r",draw=black!20,line width=0.2pt] (axis cs:"
                       + f"{phase-0.5},{row-0.5}" + r") rectangle (axis cs:"
                       + f"{phase+0.5},{row+0.5}" + ");\n")
            result += (r"\node[font=\scriptsize,text=" + ink + r"] at (axis cs:"
                       + f"{phase},{row}" + ") {" + str(value) + "};\n")
    result += r"""\draw[black,densely dashed,line width=0.6pt]
(axis cs:14.5,-0.5) -- (axis cs:14.5,18.5);
\end{axis}
\node[anchor=north,font=\scriptsize] at ([yshift=-1.13cm]panelB.south)
{19 varying coordinates; 302 constant};
\node[anchor=north,font=\scriptsize] at ([yshift=-1.52cm]panelB.south)
{$x(15)=x(11)$; no fixed-target certificate};
"""
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "larger_networks/results")
    args = parser.parse_args()
    results = args.results.resolve()
    protocol_path = ROOT / "larger_networks/protocol/contexts.json"
    protocol = json.loads(protocol_path.read_text())
    expected_ids = [spec[0] for spec in PANELS]
    if not protocol.get("frozen") or protocol["figure_context_ids"] != expected_ids:
        raise ValueError("Figure selection differs from the frozen protocol")
    contexts = {context["id"]: context for context in protocol["contexts"]}
    figures = ROOT / "figures"
    data_dir = figures / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    qa = ROOT / "qa"
    qa.mkdir(exist_ok=True)
    panels = [load_panel(results, *spec, data_dir) for spec in PANELS]
    for panel in panels:
        meta = panel["metadata"]
        if meta["input_assignment"] != contexts[meta["context_id"]]["input_assignment"]:
            raise ValueError("Plotted input assignment differs from the frozen protocol")
    if panels[1].get("kind") != "cycle":
        raise ValueError("The preselected macrophage cycle panel is required")
    source = r"""\documentclass[tikz,border=4pt]{standalone}
\usepackage[T1]{fontenc}
\usepackage{lmodern,amsmath,pgfplots}
\pgfplotsset{compat=1.18}
\definecolor{curveblue}{RGB}{26,90,137}
\definecolor{curveorange}{RGB}{190,105,28}
\pdfinfo{/Author () /Title (Support reduction and verified cycle in larger Boolean networks)}
\begin{document}
\begin{tikzpicture}
\pgfplotsset{every axis/.append style={
title style={font=\small\bfseries},
axis line style={black!50},tick label style={font=\footnotesize},
label style={font=\footnotesize}}}
PANELS
\end{tikzpicture}
\end{document}
""".replace("PANELS", support_panel_tex(panels[0]) + cycle_panel_tex(panels[1]))
    tex_path = figures / "fig5_larger_network_support.tex"
    tex_path.write_text(source)
    run = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                          tex_path.name], cwd=figures, capture_output=True, text=True)
    (qa / "figure5_build.txt").write_text(run.stdout + run.stderr)
    if run.returncode:
        raise RuntimeError(run.stdout[-6000:])
    log = tex_path.with_suffix(".log").read_text()
    metadata = {
        "status": "PASS", "panels_prespecified_before_results": True,
        "prespecified_context_ids": expected_ids,
        "protocol_sha256": digest(protocol_path),
        "selection_basis": "Scientific input identities frozen before outcome execution; no ranking or replacement based on favorable outcomes.",
        "description": "Panel B displays the independently verified native-source cycle; panel A displays the support reduction.",
        "panels": [panel["metadata"] for panel in panels],
        "tex_sha256": digest(tex_path),
        "pdf_sha256": digest(tex_path.with_suffix(".pdf")),
        "warnings": [line for line in log.splitlines()
                     if "Warning" in line or "Overfull" in line or "Underfull" in line],
    }
    (qa / "figure5_build.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"status": metadata["status"],
                      "contexts": metadata["prespecified_context_ids"],
                      "warnings": metadata["warnings"]}, indent=2))


if __name__ == "__main__":
    main()
