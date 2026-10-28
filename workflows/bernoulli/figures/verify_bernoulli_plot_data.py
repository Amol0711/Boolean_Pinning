#!/usr/bin/env python3
"""Check figure exports against independently produced statistics without importing the builder."""
from pathlib import Path
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figures"
DATA = FIG / "data"
STATS = ROOT / "results/bernoulli_statistics"


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a, b):
    assert math.isclose(float(a), float(b), rel_tol=2e-13, abs_tol=2e-14), (a, b)


def main():
    cell = {(r["mixture"], int(r["n"])): r for r in rows(STATS / "cell_summary.csv")}
    pmf = {(r["mixture"], int(r["k"])): r for r in rows(STATS / "pmf_intervals.csv")}
    surv = {(r["mixture"], int(r["threshold"])): r for r in rows(STATS / "residual_survival.csv")}
    metadata = json.loads((FIG / "bernoulli_plot_metadata.json").read_text())
    checked = {"pmf_bins": 0, "certificate_cell_records": 0, "normalized_moment_cell_records": 0,
               "residual_survival_rows_including_zero": 0}
    for name, expected in metadata["input_hashes"].items():
        assert digest(STATS / name) == expected, name
    for m in "ABCDE":
        for row in rows(DATA / ("bernoulli_pmf_" + m + ".csv")):
            k = int(row["k"])
            source = pmf[m, k]
            q = float(source["poisson_probability"])
            close(row["reference"], q)
            close(row["x"], k + metadata["pmf_bin_display_offsets"][m])
            for export, source_column in (("lo", "prob_lower"), ("hi", "prob_upper"),
                                          ("outerlo", "cp_lower"), ("outerhi", "cp_upper")):
                close(row[export], float(source[source_column]) - q)
            for export, source_column in (("prob_lower", "prob_lower"), ("prob_upper", "prob_upper"),
                                          ("probability_cp_lower", "cp_lower"), ("probability_cp_upper", "cp_upper")):
                close(row[export], source[source_column])
            close(row["err"], float(row["hi"]) - float(row["lo"]))
            close(row["outerplus"], float(row["outerhi"]) - float(row["lo"]))
            close(row["outerminus"], float(row["lo"]) - float(row["outerlo"]))
            checked["pmf_bins"] += 1
        for certificate in ("g", "g_prime"):
            for row in rows(DATA / ("bernoulli_failure_" + certificate + "_" + m + ".csv")):
                n = int(row["n"])
                source = cell[m, n]
                assert int(row["count"]) == int(source[certificate + "_failures"])
                close(row["true_rate"], int(row["count"])/int(source["M"]))
                if int(row["count"]):
                    close(row["p"], row["true_rate"])
                else:
                    assert row["p"] == "nan" and float(row["cp_lo"]) == 0
                close(row["cp_lo"], source[certificate + "_cp_lo"])
                close(row["cp_hi"], source[certificate + "_cp_hi"])
                close(row["minus"], float(row["true_rate"]) - float(row["cp_lo"]))
                close(row["plus"], float(row["cp_hi"]) - float(row["true_rate"]))
                close(row["x"], n*metadata["network_size_display_offsets"][m])
                checked["certificate_cell_records"] += 1
        for row in rows(DATA / ("bernoulli_moments_" + m + ".csv")):
            n = int(row["n"])
            source = cell[m, n]
            close(row["x"], n*metadata["network_size_display_offsets"][m])
            for j, stem, norm_key in ((1, "mean", "mu"), (2, "second_moment", "poisson_second_moment")):
                norm = float(source[norm_key])
                for prefix, suffix in (("lo", "_lower"), ("hi", "_upper"),
                                       ("outerlo", "_sampling_lo"), ("outerhi", "_sampling_hi")):
                    close(float(row[prefix+str(j)])*norm, source[stem+suffix])
                close(row["err"+str(j)], float(row["hi"+str(j)]) - float(row["lo"+str(j)]))
                close(row["outerplus"+str(j)], float(row["outerhi"+str(j)]) - float(row["lo"+str(j)]))
                close(row["outerminus"+str(j)], float(row["lo"+str(j)]) - float(row["outerlo"+str(j)]))
            checked["normalized_moment_cell_records"] += 1
        for row in rows(DATA / ("bernoulli_survival_" + m + ".csv")):
            threshold = int(row["x"])
            source = surv[m, threshold]
            close(row["survival"], source["survival"])
            if float(row["survival"]):
                close(row["plot_survival"], row["survival"])
            else:
                assert row["plot_survival"] == "nan"
            checked["residual_survival_rows_including_zero"] += 1
    assert checked["pmf_bins"] == 45
    assert checked["certificate_cell_records"] == 50
    assert checked["normalized_moment_cell_records"] == 25
    assert checked["residual_survival_rows_including_zero"] == len(surv)
    result = dict(status="PASS", comparisons=checked,
                  input_hashes_verified=True,
                  unresolved_draws_imputed=0,
                  source_metadata_sha256=digest(FIG / "bernoulli_plot_metadata.json"))
    qa = ROOT / "qa/bernoulli_figures"
    qa.mkdir(parents=True, exist_ok=True)
    (qa / "plot_data_validation.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
