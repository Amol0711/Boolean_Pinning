#!/usr/bin/env python3
"""Verify the grouped orbit display against all retained cells and coordinates.

This checks a derived display of retained records. It does not rerun the
biological model. The original full-coordinate asset and data remain intact.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

FIGURES = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def verify():
    data = FIGURES / "data"
    cycle_path = data / "fig5_mac_bacteria_cycle.csv"
    orbit_path = data / "fig5_mac_bacteria_full_orbit.csv"
    cd4_path = data / "fig5_cd4_apc.csv"
    original_path = FIGURES / "fig5_larger_network_support.tex"
    grouped_path = FIGURES / "fig5_larger_network_support_grouped.tex"
    map_path = FIGURES / "fig5_grouped_profile_map.json"
    cycles = read_csv(cycle_path)
    orbit = read_csv(orbit_path)
    cd4 = read_csv(cd4_path)
    original = original_path.read_text()
    grouped = grouped_path.read_text()
    mapping = json.loads(map_path.read_text())
    updates = [11, 12, 13, 14, 15]
    require(mapping["source"] == "data/fig5_mac_bacteria_cycle.csv", "Incorrect mapping source")
    require(mapping["source_sha256"] == digest(cycle_path), "Mapping source hash differs")
    require(mapping["updates"] == updates, "Displayed update sequence differs")
    require(len(cycles) == 19, "Expected 19 varying-coordinate records")
    require(mapping["varying_coordinates"] == 19, "Mapping coordinate total differs")

    # Derive profiles from the retained CSV in first-occurrence order.
    derived = {}
    for row in cycles:
        values = tuple(int(row[f"update_{t}"]) for t in updates)
        require(all(value in (0, 1) for value in values), "A retained state is non-Boolean")
        derived.setdefault(values, []).append({
            "node_index_zero_based": int(row["node_index_zero_based"]),
            "source_node_name": row["source_node_name"],
            "display_label": row["display_label"],
        })
    expected = [
        {"row": i, "multiplicity": len(members), "values": list(values), "members": members}
        for i, (values, members) in enumerate(derived.items())
    ]
    require(mapping["groups"] == expected, "Profile membership or multiplicity differs")
    require(len(expected) == 4, "Expected four distinct profiles")
    require([g["multiplicity"] for g in expected] == [12, 3, 3, 1], "Unexpected profile counts")
    require(sum(g["multiplicity"] for g in expected) == 19, "Profile counts do not total 19")

    # Read every displayed rectangle and its explicit Boolean digit from TikZ.
    b = grouped.split(r"\begin{axis}[name=panelB", 1)[1].split(r"\end{axis}", 1)[0]
    cells = re.findall(
        r"\\path\[fill=(curveblue|white),draw=black!20,line width=0\.2pt\] "
        r"\(axis cs:([-\d.]+),([-\d.]+)\) rectangle "
        r"\(axis cs:([-\d.]+),([-\d.]+)\);\s*"
        r"\\node\[font=\\fontsize\{8\}\{9\.5\}\\selectfont,text=(white|black!75)\] "
        r"at \(axis cs:(\d+),(\d+)\) \{([01])\};", b
    )
    require(len(cells) == 20, "The grouped raster must contain 20 explicit cells")
    rendered = {}
    for fill, left, low, right, high, ink, column, row, digit in cells:
        t, i, value = int(column), int(row), int(digit)
        require((t, i) not in rendered, "A grouped raster cell is duplicated")
        require((float(left), float(low), float(right), float(high)) ==
                (t - 0.5, i - 0.5, t + 0.5, i + 0.5), "Incorrect raster-cell bounds")
        require(fill == ("curveblue" if value else "white"), "Raster fill disagrees with digit")
        require(ink == ("white" if value else "black!75"), "Raster digit contrast is incorrect")
        rendered[t, i] = value
    require("yticklabels={12,3,3,1}" in b, "Multiplicity labels differ")
    require("ylabel={Coordinates per profile}" in b, "Multiplicity meaning is missing")
    require("xtick={11,12,13,14,15}" in b, "Raster update ticks differ")
    require("ymin=-0.5,ymax=3.5,y dir=reverse" in b, "Raster row order differs")
    require("(axis cs:14.5,-0.5) -- (axis cs:14.5,3.5)" in b,
            "First repeated-state column is not separated")
    verified_cells = 0
    for group in expected:
        for member in group["members"]:
            record = next(row for row in cycles if row["source_node_name"] == member["source_node_name"])
            for t in updates:
                require(rendered[t, group["row"]] == int(record[f"update_{t}"]),
                        "Grouped raster does not reconstruct a retained cell")
                verified_cells += 1
    require(verified_cells == 95, "Not all 95 original raster cells were checked")

    # Verify the cycle and its first repetition using the complete retained orbit.
    names = [name for name in orbit[0] if name != "synchronous_update"]
    require(len(names) == 321, "The retained orbit must have 321 coordinates")
    by_time = {int(row["synchronous_update"]): row for row in orbit}
    require(sorted(by_time) == list(range(16)), "The retained orbit update range differs")
    seen = {}
    first_repeat = None
    for t in range(16):
        state = tuple(int(by_time[t][name]) for name in names)
        if state in seen:
            first_repeat = (t, seen[state])
            break
        seen[state] = t
    require(first_repeat == (15, 11), "The first repeated state is not x(15)=x(11)")
    varying = [name for name in names if len({by_time[t][name] for t in range(11, 15)}) > 1]
    require(set(varying) == {row["source_node_name"] for row in cycles},
            "The raster does not include every varying coordinate")
    require(len(names) - len(varying) == mapping["constant_coordinates_on_cycle"] == 302,
            "Constant-coordinate scope differs")
    for row in cycles:
        name = row["source_node_name"]
        require(names[int(row["node_index_zero_based"])] == name, "A source coordinate index differs")
        for t in updates:
            require(row[f"update_{t}"] == by_time[t][name], "Raster record differs from full orbit")

    # The CD4 scientific axes and both complete plotting commands are retained.
    def axis_a(source):
        return source.split(r"\begin{axis}[name=panelA", 1)[1].split(r"\end{axis}", 1)[0]
    old_a, new_a = axis_a(original), axis_a(grouped)
    required_axis = ["xmin=-0.2,xmax=14.25,ymin=-4,ymax=168",
                     "xtick={0,2,...,14},ytick={0,40,...,160}",
                     "xlabel={Support-reduction round $k$}",
                     "ylabel={Support-bound cardinality}"]
    require(all(item in old_a and item in new_a for item in required_axis),
            "A CD4 axis range, tick sequence or axis label changed")
    plots = lambda source: re.findall(r"\\addplot.*?;", source, re.S)
    require(len(plots(old_a)) == 2 and plots(old_a) == plots(new_a),
            "A CD4 plotting command changed")
    require([int(row["support_round"]) for row in cd4] == list(range(15)),
            "CD4 support rounds differ")
    require(next(int(row["support_round"]) for row in cd4 if row["optimal_support"] == "0") == 13,
            "The controlled CD4 support does not first vanish at round 13")
    require(cd4[-1]["baseline_support"] == "17", "The uncontrolled CD4 support endpoint differs")

    return {
        "status": "PASS",
        "scope": "Derived display verification against retained records; no biological computation replay",
        "profile_multiplicities": [g["multiplicity"] for g in expected],
        "original_cells_checked": verified_cells,
        "displayed_profile_cells_checked": len(cells),
        "coordinate_names_retained_in_mapping": len(cycles),
        "orbit_coordinates": len(names),
        "constant_coordinates_on_cycle": 302,
        "first_repetition": {"update": 15, "previous_update": 11, "period": 4},
        "cd4_records": len(cd4),
        "cd4_scientific_axes_and_plot_commands_unchanged": True,
        "hashes": {str(path.relative_to(FIGURES)): digest(path) for path in
                   [cycle_path, orbit_path, cd4_path, original_path, grouped_path, map_path]},
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional verification record outside retained data")
    args = parser.parse_args()
    result = verify()
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(encoded, end="")
