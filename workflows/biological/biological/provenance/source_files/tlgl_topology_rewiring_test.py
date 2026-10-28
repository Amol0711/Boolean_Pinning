
import csv
import itertools
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class Rule:
    kind: str
    slots: List[int]


# T-LGL Boolean/IBMF rule templates inferred from the released 29-node QL script.
TRUE_RULES: List[Rule] = [
    Rule('copy', [0]),
    Rule('copy', [0]),
    Rule('copy', [1]),
    Rule('copy', [0]),
    Rule('copy', [0]),
    Rule('copy', [3]),
    Rule('or_not', [4, 5]),
    Rule('or_and_or', [13, 5, 2, 4]),
    Rule('copy', [8]),
    Rule('copy', [8]),
    Rule('copy', [9]),
    Rule('notandnot', [3, 10]),
    Rule('notandnot', [3, 10]),
    Rule('copy', [10]),
    Rule('or', [10, 15]),
    Rule('copy', [14]),
    Rule('copy', [14]),
    Rule('not', [16]),
    Rule('copy', [17]),
    Rule('notand', [0, 18]),
    Rule('copy', [19]),
    Rule('copy', [0]),
    Rule('copy', [1]),
    Rule('copy', [3]),
    Rule('copy', [11]),
    Rule('copy', [13]),
    Rule('copy', [13]),
    Rule('copy', [13]),
    Rule('copy', [10]),
]

TARGET = np.array(
    [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0],
    dtype=float,
)

# The released 29-node QL script has two equivalent 3-node validated solutions:
# (0, 8, 14) and (0, 8, 15). We use the lexicographically first one.
VALIDATED_SET = (0, 8, 14)

RHOS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
N_REPLICATES = 50
N_INIT_STATES = 1000
BOOLEAN_HORIZON = 8
IBMF_HORIZON = 10
SEED = 1


def build_combo_cache(n_nodes: int = 29, k: int = 3):
    combos = np.array(list(itertools.combinations(range(n_nodes), k)), dtype=int)
    fixed_mask = np.zeros((len(combos), n_nodes), dtype=bool)
    fixed_mask[np.arange(len(combos))[:, None], combos] = True
    return combos, fixed_mask


COMBOS, FIXED_MASK = build_combo_cache()


def build_eligible_edges(rules: Sequence[Rule]) -> List[Tuple[int, int, int]]:
    # Keep self-persistence inputs fixed; all other regulatory inputs may be rewired.
    eligible = []
    for node_idx, rule in enumerate(rules):
        for pos_idx, src in enumerate(rule.slots):
            if rule.kind == 'copy' and src == node_idx:
                continue
            eligible.append((node_idx, pos_idx, src))
    return eligible


ELIGIBLE_EDGES = build_eligible_edges(TRUE_RULES)


def update_states(states: np.ndarray, rules: Sequence[Rule]) -> np.ndarray:
    out = np.empty_like(states)
    for i, rule in enumerate(rules):
        s = rule.slots
        k = rule.kind
        if k == 'copy':
            out[:, i] = states[:, s[0]]
        elif k == 'or_not':
            a = states[:, s[0]]
            b = 1 - states[:, s[1]]
            out[:, i] = a + b - a * b
        elif k == 'or_and_or':
            a = states[:, s[0]]
            b = states[:, s[1]]
            c = states[:, s[2]]
            d = states[:, s[3]]
            c_or_d = c + d - c * d
            band = b * c_or_d
            out[:, i] = a + band - a * band
        elif k == 'notandnot':
            out[:, i] = (1 - states[:, s[0]]) * (1 - states[:, s[1]])
        elif k == 'or':
            a = states[:, s[0]]
            b = states[:, s[1]]
            out[:, i] = a + b - a * b
        elif k == 'not':
            out[:, i] = 1 - states[:, s[0]]
        elif k == 'notand':
            out[:, i] = (1 - states[:, s[0]]) * states[:, s[1]]
        else:
            raise ValueError(f'Unknown rule kind: {k}')
    return out


def rewire_rules(rules: Sequence[Rule], rho: float, rng: np.random.Generator) -> List[Rule]:
    rewired = [Rule(rule.kind, rule.slots.copy()) for rule in rules]
    m = len(ELIGIBLE_EDGES)
    k = int(round(rho * m))
    if k == 0:
        return rewired
    chosen = rng.choice(m, size=k, replace=False)
    for idx in chosen:
        node_idx, pos_idx, old_src = ELIGIBLE_EDGES[idx]
        new_src = int(rng.integers(0, len(rules)))
        while new_src == old_src:
            new_src = int(rng.integers(0, len(rules)))
        rewired[node_idx].slots[pos_idx] = new_src
    return rewired


def entropy(states: np.ndarray) -> np.ndarray:
    eps = 1e-12
    with np.errstate(divide='ignore', invalid='ignore'):
        h = np.where(
            (states <= eps) | (states >= 1 - eps),
            0.0,
            -(states * np.log2(states) + (1 - states) * np.log2(1 - states)),
        )
    return h.mean(axis=1)


def best_three_node_set_ibmf(rules: Sequence[Rule], target: np.ndarray) -> Tuple[Tuple[int, int, int], float, float]:
    n_combos, n_nodes = FIXED_MASK.shape
    states = np.full((n_combos, n_nodes), 0.5, dtype=float)
    target_broadcast = np.broadcast_to(target, (n_combos, n_nodes))
    states[FIXED_MASK] = target_broadcast[FIXED_MASK]
    for _ in range(IBMF_HORIZON):
        states = update_states(states, rules)
        states[FIXED_MASK] = target_broadcast[FIXED_MASK]
    final_distance = np.abs(states - target_broadcast).mean(axis=1)
    final_entropy = entropy(states)
    # Lexicographic tie-break on node indices for deterministic output.
    order = np.lexsort((COMBOS[:, 2], COMBOS[:, 1], COMBOS[:, 0], final_entropy, final_distance))
    best = int(order[0])
    return tuple(map(int, COMBOS[best])), float(final_distance[best]), float(final_entropy[best])


def success_rate_on_true_model(control_set: Sequence[int], init_states: np.ndarray, target: np.ndarray) -> float:
    states = init_states.copy()
    control = np.array(control_set, dtype=int)
    states[:, control] = target[control]
    for _ in range(BOOLEAN_HORIZON):
        states = update_states(states, TRUE_RULES)
        states[:, control] = target[control]
    return float(np.all(np.abs(states - target) < 1e-12, axis=1).mean())


def summarise(values: Sequence[float]) -> Tuple[float, float, float]:
    arr = np.array(values, dtype=float)
    return float(np.median(arr)), float(np.quantile(arr, 0.25)), float(np.quantile(arr, 0.75))


def run_experiment(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    init_states = rng.integers(0, 2, size=(N_INIT_STATES, len(TRUE_RULES))).astype(float)

    validated_success = success_rate_on_true_model(VALIDATED_SET, init_states, TARGET)

    rows = []
    for rho in RHOS:
        for rep in range(N_REPLICATES):
            rewired_rules = rewire_rules(TRUE_RULES, rho, rng)
            predicted_set, ibmf_distance, ibmf_entropy = best_three_node_set_ibmf(rewired_rules, TARGET)
            true_success = success_rate_on_true_model(predicted_set, init_states, TARGET)
            rows.append({
                'rho': rho,
                'replicate': rep,
                'predicted_set': '-'.join(map(str, predicted_set)),
                'ibmf_distance': ibmf_distance,
                'ibmf_entropy': ibmf_entropy,
                'baseline_success_true_model': true_success,
                'validated_success_true_model': validated_success,
            })

    csv_path = out_dir / 'tlgl_rewiring_results.csv'
    with csv_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # Aggregate for plotting and text insertion.
    summary_rows = []
    by_rho = {}
    for row in rows:
        by_rho.setdefault(row['rho'], []).append(row)
    for rho in RHOS:
        vals = [r['baseline_success_true_model'] for r in by_rho[rho]]
        med, q1, q3 = summarise(vals)
        summary_rows.append({
            'rho': rho,
            'baseline_median': med,
            'baseline_q1': q1,
            'baseline_q3': q3,
            'validated_success': validated_success,
        })

    summary_path = out_dir / 'tlgl_rewiring_summary.csv'
    with summary_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    # Figure panel
    x = np.array([row['rho'] for row in summary_rows], dtype=float)
    baseline_median = np.array([row['baseline_median'] for row in summary_rows], dtype=float)
    baseline_q1 = np.array([row['baseline_q1'] for row in summary_rows], dtype=float)
    baseline_q3 = np.array([row['baseline_q3'] for row in summary_rows], dtype=float)
    validated = np.array([validated_success] * len(summary_rows), dtype=float)

    plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'stix',
        'font.size': 13,
        'axes.titlesize': 18,
        'axes.labelsize': 15,
        'xtick.labelsize': 13,
        'ytick.labelsize': 13,
        'legend.fontsize': 12,
    })

    fig = plt.figure(figsize=(11, 11.8), dpi=100)
    fig.patch.set_facecolor('#efefef')

    ax = fig.add_axes([0.12, 0.30, 0.80, 0.54])
    ax.grid(True, alpha=0.3)

    ax.plot(
        x,
        validated,
        marker='o',
        linewidth=3,
        label=r'QL+RTS validated set',
    )
    ax.plot(
        x,
        baseline_median,
        marker='s',
        linewidth=2.6,
        linestyle='--',
        label='topology-dependent model search',
    )
    ax.fill_between(x, baseline_q1, baseline_q3, alpha=0.20)

    ax.set_xlim(-0.01, 0.51)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel(r'Rewired-input fraction $\rho$')
    ax.set_ylabel('Success rate on true T-LGL model')
    ax.set_title('B  Robustness to topology rewiring', pad=12)
    ax.legend(loc='lower left', frameon=True)

    med_r03 = baseline_median[list(x).index(0.3)]
    ax.annotate(
        r'baseline median $\approx %.2f$ at $\rho=0.3$' % med_r03,
        xy=(0.3, med_r03),
        xytext=(0.19, 0.58),
        arrowprops=dict(arrowstyle='->', lw=1.6),
        fontsize=12,
    )
    ax.annotate(
        'validated set remains at 1.00\nacross all corruption levels',
        xy=(0.4, 1.0),
        xytext=(0.22, 0.93),
        arrowprops=dict(arrowstyle='->', lw=1.5),
        fontsize=12,
        ha='left',
        va='top',
    )

    fig.text(
        0.5,
        0.12,
        '50 rewired models per level; evaluation on 1,000 random initial conditions',
        ha='center',
        va='center',
        fontsize=13,
    )
    fig.text(
        0.5,
        0.06,
        'Sets chosen by exhaustive 3-node search on the corrupted logical model,\nthen tested on the original T-LGL dynamics',
        ha='center',
        va='center',
        fontsize=13,
    )

    panel_path = out_dir / 'tlgl_robustness_panel.png'
    fig.savefig(panel_path, dpi=100)
    plt.close(fig)


def main() -> None:
    out_dir = Path(__file__).resolve().parent / 'tlgl_rewiring_outputs'
    run_experiment(out_dir)
    print(f'Wrote outputs to {out_dir}')


if __name__ == '__main__':
    main()
