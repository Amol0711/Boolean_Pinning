#!/usr/bin/env python3
"""Export and compile the Bernoulli synthetic-study vector figures.

The independent statistical analysis supplies all intervals. This script only
checks, rescales and plots them; it never replaces unresolved draws by estimates.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figures"
DATA = FIG / "data"
QA = ROOT / "qa/bernoulli_figures"
MIXTURES = "ABCDE"
SIZES = (100, 400, 1600, 6400, 25600)
COLORS = {"A": "0072B2", "B": "D55E00", "C": "009E73", "D": "CC79A7", "E": "806400"}
MARKS = {"A": "*", "B": "square*", "C": "triangle*", "D": "diamond*", "E": "pentagon*"}
OFFSETS = {"A": .86, "B": .93, "C": 1., "D": 1/.93, "E": 1/.86}
PMF_OFFSETS = {"A": -.18, "B": .18, "C": -.18, "D": .18, "E": 0.}
NAMES = ("fig1_cost_distribution", "fig2_certificate_failure", "fig3_residual_and_moments")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def preserved_hashes():
    return {str(p.relative_to(FIG)): sha(p) for p in
            [*FIG.glob("fig[456]*"), *DATA.glob("fig[456]*")] if p.is_file()}


def remove_obsolete_plot_exports():
    """Remove stale exports before creating the current plot inputs."""
    obsolete = [DATA / (stem + "_" + m + ".csv")
                for stem in ("pmf", "failure", "moments", "survival") for m in MIXTURES]
    obsolete.append(FIG / "plot_metadata.json")
    for path in obsolete:
        if path.is_file():
            path.unlink()


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(name, rows, fields=None):
    with (DATA / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


PREAMBLE = r"""\documentclass[tikz,border=2pt]{standalone}
\usepackage[T1]{fontenc}
\usepackage{lmodern,amsmath,amssymb,pgfplots}
\usepgfplotslibrary{groupplots}
\pgfplotsset{compat=1.18}
\pdfinfo{/Author () /Title (Complete-sample finite-size control and certification)}
\pgfplotsset{every axis/.append style={
font=\fontsize{7.5}{9}\selectfont,
tick label style={font=\fontsize{7.5}{9}\selectfont},
label style={font=\fontsize{7.5}{9}\selectfont},
title style={font=\fontsize{7.8}{9.4}\selectfont\bfseries},
axis line style={black!60},tick style={black!60},
grid=major,grid style={black!10},axis background/.style={fill=white},
major tick length=2pt,minor tick length=1.2pt,
legend style={draw=none,fill=none,font=\fontsize{7.5}{9}\selectfont,cells={anchor=west}},
scaled ticks=false,unbounded coords=jump}}
"""


def document(body):
    colors = "\n".join(r"\definecolor{mix" + m + "}{HTML}{" + c + "}" for m, c in COLORS.items())
    return PREAMBLE + colors + "\n\\begin{document}\n\\begin{tikzpicture}\n" + body + "\n\\end{tikzpicture}\n\\end{document}\n"


def range_plot(m, table, power="", legend=False):
    """Thin outer sampling interval, thick inner computational interval."""
    y = "lo" + power
    body = (r"\addplot[only marks,mark=none,forget plot,color=mix" + m + r"!60!white,"
            r"error bars/.cd,y dir=both,y explicit,error bar style={line width=.4pt},"
            r"error mark options={rotate=90,mark size=2.6pt,line width=.4pt}] "
            r"table[col sep=comma,x=x,y=" + y + ",y error plus=outerplus" + power +
            ",y error minus=outerminus" + power + "] {data/" + table + "};\n")
    body += (r"\addplot+[only marks,color=mix" + m + ",mark=" + MARKS[m] +
             r",mark size=1.45pt,mark options={draw=mix" + m + ",fill=mix" + m + "},"
             r"error bars/.cd,y dir=plus,y explicit,error bar style={line width=1.1pt},"
             r"error mark options={rotate=90,mark size=1.5pt,line width=1.1pt}] "
             r"table[col sep=comma,x=x,y=" + y + ",y error=err" + power + "] {data/" + table + "};\n")
    if legend:
        body += r"\addlegendentry{" + m + "}\n"
    return body


def check_nested(lo, hi, outerlo, outerhi):
    assert -1e-12 <= lo <= hi + 1e-12
    assert outerlo <= lo + 1e-12 and outerhi + 1e-12 >= hi
    assert math.isfinite(outerlo) and math.isfinite(outerhi)


def export_data(analysis):
    summary = read_csv(analysis / "cell_summary.csv")
    pmf = read_csv(analysis / "pmf_intervals.csv")
    survival = read_csv(analysis / "residual_survival.csv")
    cells = {(r["mixture"], int(r["n"])): r for r in summary}
    assert set(cells) == {(m, n) for m in MIXTURES for n in SIZES}
    assert len(summary) == 25 and sum(int(r["M"]) for r in summary) == 250000
    assert all(int(r["M"]) == 10000 for r in summary)
    assert all(int(r["g_prime_failures"]) <= int(r["g_failures"]) for r in summary)
    assert len(pmf) == 45
    limit = {1: [math.inf, -math.inf], 2: [math.inf, -math.inf]}
    difference_max = 0.0
    zero_cells = []
    residual_maxima = {}
    for m in MIXTURES:
        bins = {int(r["k"]): r for r in pmf if r["mixture"] == m and int(r["n"]) == 25600}
        assert set(bins) == set(range(9))
        out = []
        for k in range(9):
            r = bins[k]
            lo, hi, outerlo, outerhi = (float(r[key]) for key in ("prob_lower", "prob_upper", "cp_lower", "cp_upper"))
            check_nested(lo, hi, outerlo, outerhi)
            assert outerhi <= 1 + 1e-12
            reference = float(r["poisson_probability"])
            dl, du, cl, cu = (float(r[key]) for key in
                             ("difference_lower", "difference_upper", "difference_cp_lower", "difference_cp_upper"))
            assert all(abs(a-b) < 1e-13 for a, b in
                       zip((dl, du, cl, cu), (lo-reference, hi-reference, outerlo-reference, outerhi-reference)))
            difference_max = max(difference_max, abs(cl), abs(cu))
            out.append(dict(k=k, x=k+PMF_OFFSETS[m], lo=dl, hi=du, err=du-dl,
                            outerlo=cl, outerhi=cu, outerminus=dl-cl,
                            outerplus=cu-dl, reference=reference, prob_lower=lo, prob_upper=hi,
                            probability_cp_lower=outerlo, probability_cp_upper=outerhi))
        assert sum(r["prob_lower"] for r in out) <= 1+1e-12 <= sum(r["prob_upper"] for r in out)+1e-12
        write_csv("bernoulli_pmf_" + m + ".csv", out)
        moments = []
        for cert in ("g", "g_prime"):
            failures = []
            zeros = []
            for n in SIZES:
                r = cells[m, n]
                count = int(r[cert + "_failures"])
                p, lo, hi = (float(r[cert + suffix]) for suffix in ("_rate", "_cp_lo", "_cp_hi"))
                assert abs(p-count/10000) < 1e-14
                check_nested(p, p, lo, hi)
                x = n * OFFSETS[m]
                failures.append(dict(n=n, x=x, p=p if count else "nan", true_rate=p, count=count,
                                     minus=p-lo, plus=hi-p, cp_lo=lo, cp_hi=hi))
                if count == 0:
                    assert lo == 0
                    zeros.append(dict(n=n, x=x, upper=hi))
                    zero_cells.append(dict(mixture=m, n=n, certificate=cert, count=0, cp_upper=hi))
            write_csv("bernoulli_failure_" + cert + "_" + m + ".csv", failures)
            write_csv("bernoulli_zero_" + cert + "_" + m + ".csv", zeros, ["n", "x", "upper"])
        for n in SIZES:
            r = cells[m, n]
            row = dict(n=n, x=n*OFFSETS[m])
            for power, stem, normalizer in ((1, "mean", float(r["mu"])), (2, "second_moment", float(r["poisson_second_moment"]))):
                lo, hi, outerlo, outerhi = (float(r[stem + suffix])/normalizer for suffix in ("_lower", "_upper", "_sampling_lo", "_sampling_hi"))
                check_nested(lo, hi, outerlo, outerhi)
                assert outerlo >= 0
                row.update({"lo"+str(power): lo, "hi"+str(power): hi,
                            "err"+str(power): hi-lo, "outerlo"+str(power): outerlo,
                            "outerhi"+str(power): outerhi, "outerminus"+str(power): lo-outerlo,
                            "outerplus"+str(power): outerhi-lo})
                limit[power] = [min(limit[power][0], outerlo), max(limit[power][1], outerhi)]
            moments.append(row)
        write_csv("bernoulli_moments_" + m + ".csv", moments)
        rr = [dict(x=int(r["threshold"]), survival=float(r["survival"]),
                   plot_survival=float(r["survival"]) if float(r["survival"]) > 0 else "nan")
              for r in survival if r["mixture"] == m and int(r["n"]) == 25600]
        assert rr and all(0 <= r["survival"] <= 1 for r in rr)
        assert all(a["x"] < b["x"] and a["survival"] >= b["survival"] for a, b in zip(rr, rr[1:]))
        residual_maxima[m] = max(r["x"] for r in rr if r["survival"] > 0)
        write_csv("bernoulli_survival_" + m + ".csv", rr)
    return dict(moment_limits=limit, pmf_difference_extent=max(.005, math.ceil(difference_max*1.15/.005)*.005),
                zero_failure_cells=zero_cells, residual_maxima=residual_maxima)


def make_sources(record):
    extent = record["pmf_difference_extent"]
    step = .005 if extent <= .01 else .01
    nticks = int(math.floor(extent/step + 1e-10))
    ticks = ",".join(str(step*j) for j in range(-nticks, nticks+1))
    body = r"""\begin{groupplot}[group style={group size=3 by 1,horizontal sep=1.05cm},
width=4.25cm,height=4.55cm,xmin=-.3,xmax=8.35,
""" + "ymin=" + str(-extent) + ",ymax=" + str(extent) + ",ytick={" + ticks + "}," + r"""
xtick={0,2,4,6,8},xticklabels={$0$,$2$,$4$,$6$,$\geq8$},
xlabel={Minimum pin count},
yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=2},
legend style={at={(.98,.98)},anchor=north east,row sep=-2pt}]
"""
    panels = [("A. Degree--rule association", "AB"), (r"B. Common $\lambda_*$", "CD"), ("C. Zero output degree", "E")]
    for j, (title, mixtures) in enumerate(panels):
        body += r"\nextgroupplot[title={" + title + "},ylabel={" + (r"Empirical $-$ Poisson probability" if j == 0 else "") + "}]\n"
        body += r"\addplot[black!55,densely dashed,line width=.65pt,domain=-.2:8.2,samples=2,forget plot] {0};" + "\n"
        for m in mixtures:
            body += range_plot(m, "bernoulli_pmf_" + m + ".csv", legend=True)
    body += r"""\end{groupplot}
\node[font=\fontsize{7.5}{9}\selectfont,align=center,anchor=north]
at ([yshift=-1.02cm]group c2r1.south)
{Thick bars: computational ranges\qquad Thin bars: pointwise 95\% binomial intervals\\
Dashed zero line: Poisson reference};
"""
    (FIG / (NAMES[0] + ".tex")).write_text(document(body))

    body = r"""\begin{groupplot}[group style={group size=2 by 1,horizontal sep=.8cm},
width=6.5cm,height=5.25cm,xmode=log,ymode=log,
xmin=70,xmax=36000,ymin=.0000015,ymax=.7,
xtick={100,1600,25600},xticklabels={$100$,$1600$,$25600$},
minor x tick num=0,xlabel={Network size $n$},
ytick={.00001,.001,.1},legend style={column sep=6pt}]
"""
    zeros = record["zero_failure_cells"]
    for j, (cert, title) in enumerate((("g", r"A. Original certificate $G_n$"), ("g_prime", r"B. Weaker certificate $G'_n$"))):
        body += r"\nextgroupplot[title={" + title + "},ylabel={" + ("Failure probability" if j == 0 else "") + "}"
        if j == 0:
            body += r",legend to name=bernoulli_failure_legend,legend columns=6"
        body += "]\n"
        for m in MIXTURES:
            body += (r"\addplot+[color=mix" + m + ",mark=" + MARKS[m] +
                     r",mark size=1.5pt,mark options={draw=mix" + m + ",fill=mix" + m + r"},line width=.6pt,"
                     r"error bars/.cd,y dir=both,y explicit,error bar style={line width=.4pt},"
                     r"error mark options={rotate=90,mark size=1.8pt,line width=.4pt}] "
                     r"table[col sep=comma,x=x,y=p,y error plus=plus,y error minus=minus] {data/bernoulli_failure_" + cert + "_" + m + ".csv};\n")
            if j == 0:
                body += r"\addlegendentry{" + m + "}\n"
            if any(z["certificate"] == cert and z["mixture"] == m for z in zeros):
                body += (r"\addplot[only marks,forget plot,color=mix" + m + r",mark=triangle,mark size=2.1pt,"
                         r"mark options={rotate=180,draw=mix" + m + ",fill=white,line width=.7pt}] "
                         r"table[col sep=comma,x=x,y=upper] {data/bernoulli_zero_" + cert + "_" + m + ".csv};\n")
        body += r"\addplot[black!55,densely dashed,line width=.65pt,domain=100:25600,samples=2] {1/x};" + "\n"
        if j == 0:
            body += r"\addlegendentry{$n^{-1}$}" + "\n"
    body += r"""\end{groupplot}
\node[anchor=north] at ([xshift=.4cm,yshift=-1.0cm]group c1r1.south east)
{\pgfplotslegendfromname{bernoulli_failure_legend}};
\node[font=\fontsize{7.5}{9}\selectfont,align=center,anchor=north]
at ([xshift=.4cm,yshift=-1.52cm]group c1r1.south east)
"""
    body += r"{Bars: pointwise 95\% exact binomial intervals"
    if zeros:
        body += r"\\ $\triangledown$ at the upper endpoint: zero observed failures"
    body += "};\n"
    (FIG / (NAMES[1] + ".tex")).write_text(document(body))

    max_residual = int(math.ceil(max(record["residual_maxima"].values())/20)*20)
    body = r"""\begin{groupplot}[group style={group size=3 by 1,horizontal sep=1.15cm},
width=4.25cm,height=4.7cm,legend style={draw=none}]
\nextgroupplot[title={A. Residual tail},ymode=log,xmin=0,xmax=""" + str(max_residual) + r""",
ymin=.00007,ymax=1.2,xlabel={Residual threshold $r$},ylabel={Empirical $\Pr(|R|\geq r)$},
xtick={0,100,200},ytick={1,.01,.0001},
legend to name=bernoulli_moment_legend,legend columns=5,legend style={column sep=9pt}]
"""
    for m in MIXTURES:
        body += (r"\addplot[color=mix" + m + ",mark=" + MARKS[m] + r",mark size=.9pt,mark options={draw=mix" +
                 m + ",fill=mix" + m + r"},mark repeat=15,line width=.65pt,const plot] "
                 r"table[col sep=comma,x=x,y=plot_survival] {data/bernoulli_survival_" + m + ".csv};\n")
        body += r"\addlegendentry{" + m + "}\n"
    for power, title, label in ((1, "B. Mean cost", r"$\widehat m_1/\mu$"), (2, "C. Second moment", r"$\widehat m_2/(\mu+\mu^2)$")):
        low, high = record["moment_limits"][power]
        span = max(.15, high-low)
        ymin = max(0., min(.9, low-.1*span))
        ymax = max(1.1, high+.1*span)
        body += (r"\nextgroupplot[title={" + title + r"},xmode=log,xmin=70,xmax=37000,minor x tick num=0,"
                 r"xtick={100,1600,25600},xticklabels={$100$,$1600$,$25600$},xlabel={Network size $n$},ylabel={" +
                 label + "},ymin=" + str(ymin) + ",ymax=" + str(ymax) + "]\n")
        body += r"\addplot[black!60,densely dashed,domain=80:32000,samples=2] {1};" + "\n"
        for m in MIXTURES:
            body += range_plot(m, "bernoulli_moments_" + m + ".csv", str(power))
    body += r"""\end{groupplot}
\node[anchor=north] at ([yshift=-1.0cm]group c2r1.south) {\pgfplotslegendfromname{bernoulli_moment_legend}};
\node[font=\fontsize{7.5}{9}\selectfont,align=center,anchor=north]
at ([yshift=-1.49cm]group c2r1.south)
{Moment bars\quad Thick: computational ranges\quad Thin: approximate 95\% sampling envelopes};
"""
    (FIG / (NAMES[2] + ".tex")).write_text(document(body))


def compile_figures():
    QA.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, SOURCE_DATE_EPOCH="1792022400", FORCE_SOURCE_DATE="1")
    outputs = {}
    for name in NAMES:
        with tempfile.TemporaryDirectory(prefix="boolean_bernoulli_" + name + "_") as directory:
            temporary = Path(directory)
            shutil.copy2(FIG / (name + ".tex"), temporary)
            shutil.copytree(DATA, temporary / "data")
            for _ in range(3):
                run = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", name + ".tex"],
                                     cwd=temporary, env=env, capture_output=True, text=True)
                if run.returncode:
                    raise RuntimeError(run.stdout[-5000:] + run.stderr)
            log = (temporary / (name + ".log")).read_text()
            (QA / (name + ".log")).write_text(log)
            assert "Overfull" not in log and "undefined references" not in log, name
            shutil.copy2(temporary / (name + ".pdf"), FIG / (name + ".pdf"))
            subprocess.run(["pdftoppm", "-png", "-singlefile", "-r", "180", str(FIG / (name + ".pdf")), str(QA / name)], check=True)
            outputs[name + ".pdf"] = sha(FIG / (name + ".pdf"))
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=ROOT / "results/bernoulli_statistics")
    parser.add_argument("--compile-only", action="store_true")
    args = parser.parse_args()
    args.analysis_dir = args.analysis_dir.resolve()
    DATA.mkdir(parents=True, exist_ok=True)
    unchanged = preserved_hashes()
    if not args.compile_only:
        record = export_data(args.analysis_dir)
        make_sources(record)
        record.update({
            "source_directory": str(args.analysis_dir.relative_to(ROOT)),
            "input_hashes": {name: sha(args.analysis_dir / name) for name in
                             ("cell_summary.csv", "pmf_intervals.csv", "residual_survival.csv", "summary.json")},
            "sample_accounting": {"cells": 25, "draws_per_cell": 10000, "total_draws": 250000, "excluded_draws": 0},
            "network_size_display_offsets": OFFSETS, "pmf_bin_display_offsets": PMF_OFFSETS,
            "pmf_bins": "0 through7 and one >=8 tail bin",
            "pmf_display": "Empirical bin probability minus the specified Poisson bin probability; differences and CP outer differences are read from independent statistics and verified by subtraction.",
            "interval_semantics": {"thick": "Deterministic complete-sample computational bounds",
                                   "pmf_thin": "Pointwise exact binomial outer interval supplied by independent statistics",
                                   "moment_thin": "Approximate95% sampling envelope supplied by independent statistics",
                                   "failure": "Pointwise two-sided95% Clopper-Pearson"},
            "zero_failure_display": "No point estimate is placed on log axis. An open downward triangle marks the exact CP upper endpoint; count0 and interval lower0 remain in data.",
            "moment_y_axes": "Linear; full nonnegative outer envelopes displayed",
            "poisson_normalization": "mu=-log(1-lambda_star), second raw moment mu+mu^2",
        })
        (FIG / "bernoulli_plot_metadata.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    outputs = compile_figures()
    if not args.compile_only:
        remove_obsolete_plot_exports()
    assert unchanged == preserved_hashes()
    (QA / "build.json").write_text(json.dumps({"source_date_epoch": 1792022400, "latex_passes": 3, "outputs": outputs,
        "plot_data_hashes": {p.name: sha(p) for p in sorted(DATA.glob("bernoulli_*"))}}, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "compiled", "outputs": outputs}, indent=2))


if __name__ == "__main__":
    main()
