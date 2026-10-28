#!/usr/bin/env python3
"""Export and compile the three Baseline PGFPlots figures from complete analysis tables.

Run from any directory. Requires scipy, numpy, pdflatex and PGFPlots. Sources and
small plotting tables are regenerated; PDF builds use isolated temporary folders.
No numerical outcome is estimated from a plot and no raw draw is excluded.
"""
from __future__ import annotations
import csv
import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile

from scipy.stats import poisson

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
DATA = FIG / "data"
QA = ROOT / "qa" / "figures"
ANALYSIS = ROOT / "results"
COLORS = {"A": "0072B2", "B": "D55E00", "C": "009E73", "D": "CC79A7", "E": "806400"}
MARKS = {"A": "*", "B": "square*", "C": "triangle*", "D": "diamond*", "E": "pentagon*"}
OFFSETS = {"A": .88, "B": .94, "C": 1., "D": 1/.94, "E": 1/.88}


def read_csv(name):
    with (ANALYSIS / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_table(name, rows):
    with (DATA / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


PREAMBLE = r"""\documentclass[tikz,border=2pt]{standalone}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{amsmath}
\usepackage{pgfplots}
\usepgfplotslibrary{groupplots}
\pgfplotsset{compat=1.18}
\usepackage{hyperref}
\hypersetup{pdfauthor={},pdftitle={Baseline finite-size numerical study}}
\pgfplotsset{every axis/.append style={
  font=\small,tick label style={font=\footnotesize},
  label style={font=\small},title style={font=\small},
  axis line style={black!65},tick style={black!65},
  grid=major,grid style={black!10},axis background/.style={fill=white},
  major tick length=2.4pt,minor tick length=1.4pt,
  legend style={draw=none,fill=none,font=\footnotesize,cells={anchor=west}},
  scaled ticks=false,unbounded coords=jump
}}
"""


def document(body):
    colors = "\n".join(r"\definecolor{mix"+m+r"}{HTML}{"+c+"}" for m,c in COLORS.items())
    return PREAMBLE + colors + "\n\\begin{document}\n\\begin{tikzpicture}\n" + body + "\n\\end{tikzpicture}\n\\end{document}\n"


def empirical_plot(m, table, x="x", y="lo", err="err", extra="", legend=False):
    code = (r"\addplot+[only marks,color=mix"+m+",mark="+MARKS[m]+r",mark size=1.7pt,mark options={draw=mix"+m+",fill=mix"+m+"},"
            r"error bars/.cd,y dir=plus,y explicit,error bar style={line width=.7pt},"
            r"error mark options={rotate=90,mark size=2pt,line width=.7pt},/pgfplots/.cd,"+extra+"] "
            r"table[col sep=comma,x="+x+",y="+y+",y error="+err+"] {data/"+table+"};\n")
    if legend:
        code += r"\addlegendentry{"+m+"}\n"
    return code


def build_sources():
    for path in (FIG, DATA, QA):
        path.mkdir(parents=True,exist_ok=True)
    protocol = json.loads((ROOT / "Protocol_Baseline.json").read_text())
    summary = read_csv("cell_summary.csv")
    dist = read_csv("cost_distribution.csv")
    surv = read_csv("residual_survival.csv")
    maximum_n = max(protocol["sizes"])
    assert len(summary) == protocol["total_cells"]
    assert sum(int(r["draws"]) for r in summary) == protocol["total_networks"]
    by_cell = {(r["mixture"],int(r["n"])):r for r in summary}
    largest_R = 0
    moment_limits={1:[float("inf"),0.],2:[float("inf"),0.]}
    for m in protocol["mixtures"]:
        lam = protocol["mixtures"][m]["lambda"]
        mu = -math.log1p(-lam)
        rows = {int(r["k"]): r for r in dist if r["mixture"]==m and int(r["n"])==maximum_n}
        max_k = max(rows)
        pmf = []
        for k in range(9):
            if k < 8:
                lo = float(rows[k]["pmf_lower"]) if k in rows else 0.
                hi = float(rows[k]["pmf_upper"]) if k in rows else 0.
                reference = float(poisson.pmf(k,mu))
            else:
                r = rows[min(7,max_k)]
                lo,hi = 1-float(r["cdf_upper"]),1-float(r["cdf_lower"])
                reference = float(poisson.sf(7,mu))
            assert -1e-14 <= lo <= hi+1e-14 <= 1+1e-14
            pmf.append(dict(x=k,lo=max(0,lo),hi=max(0,hi),err=max(0,hi-lo),reference=reference))
        write_table("pmf_"+m+".csv",pmf)
        failures, moments = [], []
        for n in protocol["sizes"]:
            r=by_cell[m,n]
            p=float(r["certificate_failure_fraction"])
            cp_lo=float(r["certificate_failure_cp95_lower"])
            cp_hi=float(r["certificate_failure_cp95_upper"])
            failures.append(dict(x=n,p=p if p>0 else "nan",minus=p-cp_lo,plus=cp_hi-p,upper=cp_hi,count=r["certificate_failure_count"]))
            a=float(r["raw_moment_1_lower"])/mu
            b=float(r["raw_moment_1_upper"])/mu
            c=float(r["raw_moment_2_lower"])/(mu+mu*mu)
            d=float(r["raw_moment_2_upper"])/(mu+mu*mu)
            moment_limits[1]=[min(moment_limits[1][0],a),max(moment_limits[1][1],b)]
            moment_limits[2]=[min(moment_limits[2][0],c),max(moment_limits[2][1],d)]
            moments.append(dict(x=n*OFFSETS[m],n=n,lo1=a,hi1=b,err1=b-a,lo2=c,hi2=d,err2=d-c))
        write_table("failure_"+m+".csv",failures)
        write_table("moments_"+m+".csv",moments)
        rr=[dict(x=int(r["threshold"]),survival=float(r["residual_survival_ge"]))
            for r in surv if r["mixture"]==m and int(r["n"])==maximum_n and float(r["residual_survival_ge"])>0]
        largest_R=max(largest_R,max(r["x"] for r in rr))
        write_table("survival_"+m+".csv",rr)

    # Figure 1: exact disjoint tail bin, with interval lower endpoints marked.
    body=r"""\begin{groupplot}[group style={group size=3 by 1,horizontal sep=1.1cm},
width=5.2cm,height=5.0cm,xmin=-.3,xmax=8.35,ymin=0,ymax=.66,
xtick={0,2,4,6,8},xticklabels={$0$,$2$,$4$,$6$,$\geq8$},
ytick={0,.2,.4,.6},xlabel={Minimum pin count},ylabel={Probability},
legend style={at={(.98,.97)},anchor=north east,row sep=-2pt}]
"""
    panels=[("(a) Degree--rule association",["A","B"]),("(b) Common $\\lambda_*=0.7$",["C","D"]),("(c) Zero output degree",["E"])]
    for title,mixtures in panels:
        body+=r"\nextgroupplot[title={"+title+"}]\n"
        for m in mixtures:
            body+=r"\addplot[color=mix"+m+r",line width=.8pt,dashed,forget plot] table[col sep=comma,x=x,y=reference] {data/pmf_"+m+".csv};\n"
            body+=empirical_plot(m,"pmf_"+m+".csv")
            label=m+r" ($\lambda_*="+str(protocol["mixtures"][m]["lambda"])+"$)"
            body+=r"\addlegendentry{"+label+"}\n"
    body+=r"""\end{groupplot}
\node[font=\footnotesize,anchor=north] at ([yshift=-1.15cm]group c2r1.south)
{Markers and upward bars: empirical bounds\qquad Dashed curves: Poisson reference};
"""
    (FIG/"fig1_cost_distribution.tex").write_text(document(body))

    # Figure 2: actual nonzero proportions and exact pointwise binomial CIs.
    body=r"""\begin{axis}[width=16.2cm,height=6.4cm,xmode=log,ymode=log,
xmin=75,xmax=34500,ymin=.0000015,ymax=.7,
xtick={100,400,1600,6400,25600},xticklabels={$100$,$400$,$1600$,$6400$,$25600$},
minor x tick num=0,xlabel={Network size $n$},ylabel={Certificate failure probability},
legend style={at={(.02,.03)},anchor=south west,legend columns=3,column sep=8pt,row sep=1pt}]
"""
    for m in protocol["mixtures"]:
        failures=[by_cell[m,n] for n in protocol["sizes"]]
        body+=r"\addplot+[color=mix"+m+",mark="+MARKS[m]+r",mark size=2pt,mark options={draw=mix"+m+",fill=mix"+m+r"},line width=.8pt,error bars/.cd,y dir=both,y explicit,error bar style={line width=.65pt},error mark options={rotate=90,mark size=2pt,line width=.65pt}] table[col sep=comma,x=x,y=p,y error plus=plus,y error minus=minus] {data/failure_"+m+".csv};\n"
        body+=r"\addlegendentry{"+m+"}\n"
        for r in failures:
            if int(r["certificate_failure_count"])==0:
                n,hi=int(r["n"]),float(r["certificate_failure_cp95_upper"])
                body+=r"\addplot[only marks,forget plot,color=mix"+m+r",mark=triangle*,mark options={rotate=180}] coordinates {("+str(n)+","+str(hi)+")};\n"
                body+=r"\node[color=mix"+m+r",font=\scriptsize,anchor=south west] at (axis cs:"+str(n)+","+str(hi)+") {"+m+": 0};\n"
    body+=r"""\addplot[black!65,densely dashed,line width=.8pt,domain=100:25600,samples=2] {1/x};
\addlegendentry{$n^{-1}$ reference}
\end{axis}
"""
    (FIG/"fig2_certificate_failure.tex").write_text(document(body))

    # Figure 3: residual survival and conservative moment identification ranges.
    body=r"""\begin{groupplot}[group style={group size=3 by 1,horizontal sep=1.35cm},
width=5.1cm,height=5.3cm,ymode=log,legend style={draw=none}]
\nextgroupplot[title={(a) Residual survival},xmin=0,xmax="""+str(math.ceil(largest_R/20)*20)+r""",
ymin=.00007,ymax=1.2,xlabel={Residual threshold $r$},ylabel={Empirical $\Pr(|R|\geq r)$},
xtick={0,50,100,150,200},ytick={1,.01,.0001},
legend to name=baseline_shared_legend,legend columns=5,
legend style={column sep=8pt}]
"""
    for m in protocol["mixtures"]:
        body+=r"\addplot[color=mix"+m+",mark="+MARKS[m]+r",mark size=1.1pt,mark options={draw=mix"+m+",fill=mix"+m+r"},mark repeat=15,line width=.7pt,const plot] table[col sep=comma,x=x,y=survival] {data/survival_"+m+".csv};\n"
        body+=r"\addlegendentry{"+m+"}\n"
    for power,title,label in [(1,"(b) First cost moment",r"$\widehat m_1/\mu$"),(2,"(c) Second cost moment",r"$\widehat m_2/(\mu+\mu^2)$")]:
        ymin=min(.8,moment_limits[power][0]*.85)
        ymax=max(1.2,moment_limits[power][1]*1.2)
        tick_values=[v for v in [.25,.5,1,2,5,10,20,50,100,200] if ymin<=v<=ymax]
        ticks="{"+",".join(str(v) for v in tick_values)+"}"
        ticklabels="{"+",".join("$"+str(v)+"$" for v in tick_values)+"}"
        body+=r"\nextgroupplot[title={"+title+r"},xmode=log,xmin=70,xmax=37000,minor x tick num=0,"
        body+=r"xtick={100,1600,25600},xticklabels={$100$,$1600$,$25600$},xlabel={Network size $n$},ylabel={"+label+"},ymin="+str(ymin)+",ymax="+str(ymax)+",ytick="+ticks+",yticklabels="+ticklabels+"]\n"
        body+=r"\addplot[black!70,densely dashed,domain=80:32000,samples=2] {1};"+"\n"
        for m in protocol["mixtures"]:
            body+=empirical_plot(m,"moments_"+m+".csv",y="lo"+str(power),err="err"+str(power))
    body+=r"""\end{groupplot}
\node[anchor=north] at ([yshift=-1.1cm]group c2r1.south) {\pgfplotslegendfromname{baseline_shared_legend}};
"""
    (FIG/"fig3_residual_and_moments.tex").write_text(document(body))
    metadata=dict(data_source=str(ANALYSIS.relative_to(ROOT)),input_hashes={name:hashlib.sha256((ANALYSIS/name).read_bytes()).hexdigest() for name in ["cell_summary.csv","cost_distribution.csv","residual_survival.csv"]},largest_n=maximum_n,
        pmf_bins="0 through 7, then a single >=8 tail bin",moment_normalization="mu=-log(1-lambda); first divided by mu, second by mu+mu^2",
        moment_x_offsets=OFFSETS,intervals="computational identification bounds, not confidence intervals",
        failure_intervals="pointwise two-sided 95% Clopper-Pearson",excluded_draws=0,
        cost_bound_stage="deterministic post-run refinement" if ANALYSIS.name=="refined" else "original predeclared resource-capped bounds")
    (FIG/"plot_metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")


def compile_figures():
    for source in sorted(FIG.glob("fig[123]_*.tex")):
        with tempfile.TemporaryDirectory(prefix="boolean_"+source.stem+"_") as directory:
            scratch=Path(directory)
            for run in range(3):
                command=["pdflatex","-interaction=nonstopmode","-halt-on-error","-output-directory",str(scratch),source.name]
                result=subprocess.run(command,cwd=FIG,capture_output=True,text=True)
                (QA/(source.stem+f"_pass{run+1}.txt")).write_text(result.stdout+result.stderr)
                if result.returncode:
                    raise RuntimeError(source.name+" compilation failed; see qa/figures")
            shutil.copy2(scratch/(source.stem+".pdf"),source.with_suffix(".pdf"))
            shutil.copy2(scratch/(source.stem+".log"),QA/(source.stem+".log"))
        subprocess.run(["pdftoppm","-png","-r","170","-singlefile",str(source.with_suffix(".pdf")),str(QA/source.stem)],check=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir",type=Path,default=ROOT/"results"/"refined",help="Directory containing complete analyzed CSV tables; default results/refined")
    args=parser.parse_args()
    ANALYSIS=args.analysis_dir.resolve()
    build_sources()
    compile_figures()
    print("Exported real-data tables and compiled all three figures.")
