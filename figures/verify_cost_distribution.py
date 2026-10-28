#!/usr/bin/env python3
"""Verify retained probability tables and their plotted coordinates.

Requires Python 3.10 or newer and SciPy (python -m pip install scipy).
This recalculates aggregate Clopper--Pearson endpoints and specified Poisson
probabilities. It does not reconstruct individual networks or replay the
experiment. The optional baseline argument adds byte-preservation comparisons.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re

try:
    from scipy.stats import beta, poisson
except ImportError as exc:
    raise SystemExit("SciPy is required for this check. Install it with: python -m pip install scipy") from exc

M = 10000
TYPES = {
    'A': ((1, .50), (3, .10)), 'B': ((1, .10), (3, .50)),
    'C': ((1, .80), (3, .20)), 'D': ((1, .35), (3, .35)),
    'E': ((0, .80), (4, .30)),
}
OFFSETS = {'A': -.18, 'B': .18, 'C': -.18, 'D': .18, 'E': 0.}
MARKS = {'A': '*', 'B': 'square*', 'C': 'triangle*', 'D': 'diamond*', 'E': 'pentagon*'}
# Recorded complete-sample bin counts, not generated observations.
LOWER_COUNTS = {
    'A': [5974, 3084, 787, 134, 18, 1, 2, 0, 0],
    'B': [2061, 3213, 2594, 1353, 554, 157, 51, 13, 0],
    'C': [2965, 3684, 2143, 879, 247, 59, 20, 3, 0],
    'D': [3047, 3611, 2166, 850, 254, 54, 15, 3, 0],
    'E': [3929, 3728, 1681, 521, 110, 28, 2, 1, 0],
}
UPPER_COUNTS = dict(LOWER_COUNTS, B=[2061, 3216, 2598, 1355, 555, 158, 51, 13, 0])

def require(condition, detail):
    if not condition:
        raise AssertionError(detail)

def close(a, b, label, tolerance=2e-12):
    require(math.isfinite(float(a)) and math.isfinite(float(b)), label + ': nonfinite')
    require(abs(float(a)-float(b)) <= tolerance, f'{label}: {a} != {b}')

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def rows(path):
    with path.open(newline='') as stream:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]

def scalar_options(text):
    return {key: float(value) for key, value in re.findall(
        r'\b(xmin|xmax|ymin|ymax)\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)', text)}

def inspect_plot(path, data, view):
    text = path.read_text()
    chunks = text.split('\\nextgroupplot')
    require(len(chunks) == 4, 'The main cost figure must have exactly three group panels')
    require('ymode=log' not in text.replace(' ', ''), 'Probability/difference panels must use linear ordinates')
    global_bounds = scalar_options(chunks[0])
    result = []
    pattern = re.compile(
        r'\\addplot\+?\s*\[([^\]]*)\]\s*table\s*\[([^\]]*)\]\s*\{([^}]+)\}\s*;', re.S)
    for index, (chunk, mixtures) in enumerate(zip(chunks[1:], ('AB', 'CD', 'E'))):
        panel_options = chunk.split('\\addplot', 1)[0]
        bounds = dict(global_bounds, **scalar_options(panel_options))
        require(set(bounds) >= {'xmin', 'xmax', 'ymin', 'ymax'}, 'Explicit plot limits required')
        absolute = view == 'mixed' and index == 0
        ordinate = 'prob_lower' if absolute else 'lo'
        found = {(m, kind): 0 for m in mixtures for kind in ('computational', 'sampling')}
        theories = {m: 0 for m in mixtures}
        all_coordinates = []
        commands = []
        for options, table, name in pattern.findall(chunk):
            match = re.fullmatch(r'data/bernoulli_pmf_([A-E])\.csv', name)
            require(match is not None, 'Unexpected table in probability panel: ' + name)
            mixture = match.group(1)
            require(mixture in mixtures, 'Mixture assigned to the wrong panel')
            columns = dict(re.findall(r'(x|y|y error plus|y error minus|y error)\s*=\s*([^,\]]+)', table))
            columns = {k: v.strip() for k, v in columns.items()}
            compact_options = re.sub(r'\s+', '', options)
            require('restrict' not in table + options and 'filter' not in table + options,
                    'No bin filtering is allowed')
            if columns.get('y') == 'reference':
                require(absolute and columns.get('x') == 'x',
                        'Poisson masses must align with the disclosed mixture display offsets')
                require('color=mix'+mixture in compact_options, 'Changed reference color identity')
                require(('solid' if mixture == 'A' else 'denselydashed') in compact_options,
                        'Changed reference line identity')
                theories[mixture] += 1
                all_coordinates.extend((r['x'], r['reference']) for r in data[mixture])
                commands.append({'mixture': mixture, 'kind': 'poisson', 'columns': columns})
                continue
            require(columns.get('x') == 'x' and columns.get('y') == ordinate,
                    f'Wrong empirical coordinates in panel {index+1}: {columns}')
            require('color=mix'+mixture in compact_options, 'Changed mixture color identity')
            if columns.get('y error') == 'err':
                kind = 'computational'
                require('ydir=plus' in compact_options, 'Computational bar must start at its lower endpoint')
                require('mark='+MARKS[mixture]+',' in compact_options,
                        'Changed mixture marker identity')
            else:
                kind = 'sampling'
                require(columns.get('y error plus') == 'outerplus' and
                        columns.get('y error minus') == 'outerminus', 'Wrong CP outer-error columns')
                require('ydir=both' in compact_options, 'CP envelope requires both directions')
            found[mixture, kind] += 1
            for row in data[mixture]:
                base = row[ordinate]
                target_low = row['probability_cp_lower'] if absolute else row['outerlo']
                target_high = row['probability_cp_upper'] if absolute else row['outerhi']
                if kind == 'computational':
                    upper = row['prob_upper'] if absolute else row['hi']
                    close(base+row['err'], upper, 'Plotted computational upper endpoint')
                    all_coordinates += [(row['x'], base), (row['x'], upper)]
                else:
                    close(base-row['outerminus'], target_low, 'Plotted CP lower endpoint')
                    close(base+row['outerplus'], target_high, 'Plotted CP upper endpoint')
                    all_coordinates += [(row['x'], target_low), (row['x'], target_high)]
            commands.append({'mixture': mixture, 'kind': kind, 'columns': columns})
        require(all(count == 1 for count in found.values()), 'Each mixture needs one complete inner and outer range')
        require(all(count == (1 if absolute else 0) for count in theories.values()),
                'Wrong number of specified Poisson curves')
        require(all(bounds['xmin'] <= x <= bounds['xmax'] and
                    bounds['ymin']-1e-12 <= y <= bounds['ymax']+1e-12 for x, y in all_coordinates),
                'A bin or uncertainty endpoint is clipped by the axes')
        if absolute:
            close(bounds['ymin'], 0., 'Absolute probability lower axis')
            close(bounds['ymax'], .65, 'Absolute probability upper axis')
        if not absolute:
            close(bounds['ymin'], -.02, 'Preserved difference lower axis')
            close(bounds['ymax'], .02, 'Preserved difference upper axis')
        require('\\geq8' in chunks[0] or '\\geq8' in chunk, 'The final cost bin must be labeled >=8')
        result.append({'panel': index+1, 'mixtures': list(mixtures), 'ordinate': ordinate,
                       'bounds': bounds, 'table_mappings': commands,
                       'all_nine_bins_and_all_endpoints_visible': True})
    return result



def check(figures, baseline=None):
    presentation_path = figures / 'cost_distribution_presentation.json'
    metadata = json.loads(presentation_path.read_text())
    require(metadata['network_size'] == 25600 and metadata['sample_size_per_mixture'] == M,
            'Changed size or complete-sample denominator')
    require(metadata['display_offsets'] == OFFSETS, 'Metadata offsets disagree')
    require(metadata['bins'] == '0 through 7 and the complete >=8 tail', 'Changed disjoint-bin semantics')
    require(metadata['empirical_abscissa'] == 'x' and metadata['diagnostic'] == 'fig1_cost_differences',
            'Metadata uses unexpected plot coordinates or diagnostic')
    data, hashes, mixture_results, outside = {}, {}, {}, []
    maximum_cp_error = maximum_poisson_error = 0.
    for mixture, types in TYPES.items():
        path = figures / 'data' / ('bernoulli_pmf_'+mixture+'.csv')
        series = data[mixture] = rows(path)
        hashes[path.name] = digest(path)
        require([row['k'] for row in series] == list(range(9)), 'Exactly nine disjoint bins are required')
        lam = sum(degree*q for degree, q in types)/2
        mu = -math.log1p(-lam)
        masses = [poisson.pmf(k, mu) for k in range(8)] + [poisson.sf(7, mu)]
        close(sum(masses), 1., 'Sum of specified Poisson bin probabilities')
        lower_counts, upper_counts = [], []
        for k, (row, mass) in enumerate(zip(series, masses)):
            lower, upper = row['prob_lower'], row['prob_upper']
            lo_count, hi_count = round(M*lower), round(M*upper)
            close(M*lower, lo_count, 'Integer definite count', 1e-8)
            close(M*upper, hi_count, 'Integer possible count', 1e-8)
            require(0 <= lo_count <= hi_count <= M, 'Invalid bin counts')
            lower_counts.append(lo_count); upper_counts.append(hi_count)
            cp_lo = 0. if lo_count == 0 else beta.ppf(.025, lo_count, M-lo_count+1)
            cp_hi = 1. if hi_count == M else beta.ppf(.975, hi_count+1, M-hi_count)
            maximum_cp_error = max(maximum_cp_error, abs(cp_lo-row['probability_cp_lower']),
                                   abs(cp_hi-row['probability_cp_upper']))
            maximum_poisson_error = max(maximum_poisson_error, abs(mass-row['reference']))
            close(cp_lo, row['probability_cp_lower'], 'Recomputed CP lower endpoint')
            close(cp_hi, row['probability_cp_upper'], 'Recomputed CP upper endpoint')
            close(mass, row['reference'], 'Law-derived Poisson probability')
            require(cp_lo <= lower+1e-12 <= upper+1e-12 <= cp_hi+1e-12, 'Non-nested ranges')
            close(row['x'], k+OFFSETS[mixture], 'Visibility-only bin offset')
            for translated, original in [('lo','prob_lower'), ('hi','prob_upper'),
                                         ('outerlo','probability_cp_lower'), ('outerhi','probability_cp_upper')]:
                close(row[translated], row[original]-row['reference'], 'Probability translation')
            close(row['err'], upper-lower, 'Computational interval length')
            close(row['outerminus'], lower-row['probability_cp_lower'], 'Lower CP distance')
            close(row['outerplus'], row['probability_cp_upper']-lower, 'Upper CP distance from lower marker')
            if not cp_lo <= mass <= cp_hi:
                outside.append({'mixture':mixture, 'bin':k, 'poisson_probability':mass,
                                'cp_lower':cp_lo, 'cp_upper':cp_hi})
        require(lower_counts == LOWER_COUNTS[mixture] and upper_counts == UPPER_COUNTS[mixture],
                'The recorded aggregate counts changed')
        require(sum(lower_counts) <= M <= sum(upper_counts), 'Invalid complete-sample bracketing')
        mixture_results[mixture] = {
            'lambda_star':lam, 'mu':mu, 'definite_counts':lower_counts, 'possible_counts':upper_counts,
            'definite_count_sum':sum(lower_counts), 'possible_count_sum':sum(upper_counts),
            'zero_cost_probability':series[0]['prob_lower'],
            'zero_cost_cp':[series[0]['probability_cp_lower'],series[0]['probability_cp_upper']],
            'zero_cost_poisson_reference':masses[0], 'poisson_tail_ge8':masses[8]}
    require(metadata['plot_data_sha256'] == hashes, 'Presentation metadata/data hashes disagree')
    for panel, mixtures, ordinate, limits in [('a','AB','prob_lower',[0,.65]),
                                            ('b','CD','lo',[-.02,.02]), ('c','E','lo',[-.02,.02])]:
        record = metadata['panels'][panel]
        require(record['mixtures'] == mixtures and record['empirical_ordinate'] == ordinate
                and record['y_limits'] == limits, 'Presentation panel metadata disagrees')
    require(metadata['panels']['a']['reference_abscissa'] == 'x'
            and metadata['panels']['a']['reference_ordinate'] == 'reference', 'Incorrect Poisson reference mapping')
    require(metadata['panels']['b']['reference'] == metadata['panels']['c']['reference'] == 0,
            'Difference reference must be zero')
    mappings = inspect_plot(figures/'fig1_cost_distribution.tex', data, 'mixed')
    diagnostic_mappings = inspect_plot(figures/'fig1_cost_differences.tex', data, 'difference')
    require(len(outside) == 1 and outside[0]['mixture'] == 'A' and outside[0]['bin'] == 6,
            'The retained reference/envelope comparison changed')
    baseline_checks = []
    if baseline is not None:
        pairs = [(figures/'data'/name, baseline/'data'/name) for name in hashes]
        pairs += [(figures/('fig1_cost_differences.'+suffix),
                   baseline/('fig1_cost_distribution.'+suffix)) for suffix in ('tex','pdf')]
        for current, original in pairs:
            require(current.read_bytes() == original.read_bytes(), 'Byte preservation failed: '+current.name)
            baseline_checks.append(current.name)
    return {
        'status':'PASS',
        'scope':'Retained aggregate endpoint and plot-source verification; no raw-network replay',
        'network_size':25600, 'draws_per_mixture':M, 'bins_per_mixture':9, 'total_bins_checked':45,
        'maximum_cp_endpoint_error':maximum_cp_error,
        'maximum_poisson_probability_error':maximum_poisson_error,
        'data_hashes_sha256':hashes, 'presentation_metadata_sha256':digest(presentation_path),
        'plot_source_sha256':digest(figures/'fig1_cost_distribution.tex'),
        'diagnostic_source_sha256':digest(figures/'fig1_cost_differences.tex'),
        'mixtures':mixture_results, 'main_plot_mappings':mappings,
        'diagnostic_plot_mappings':diagnostic_mappings,
        'poisson_references_outside_pointwise_cp_envelopes':outside,
        'baseline_byte_preservation_checked':baseline is not None,
        'baseline_files_checked':baseline_checks,
        'limitations':[
            'Pointwise finite-size probability envelopes do not provide simultaneous coverage or Poisson-approximation error bounds.',
            'Aggregate tables do not identify individual networks, cost intervals, certificates or witnesses.',
            'Text and rendered-figure legibility require visual inspection.']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-figures', type=Path,
                        help='Optional original figures directory for CSV and preserved-diagnostic byte checks')
    parser.add_argument('--output', type=Path, help='Result JSON destination; defaults to source-root qa/')
    arguments = parser.parse_args()
    figures = Path(__file__).resolve().parent
    output = arguments.output or figures.parent/'qa/cost_distribution_verification.json'
    result = check(figures, arguments.baseline_figures)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ['status', 'total_bins_checked',
          'maximum_cp_endpoint_error', 'maximum_poisson_probability_error',
          'baseline_byte_preservation_checked']}, indent=2))
