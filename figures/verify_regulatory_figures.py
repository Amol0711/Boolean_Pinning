#!/usr/bin/env python3
"""Verify retained regulatory data and the native figure mappings (stdlib only).

This is a retained-data and figure-consistency check. It does not replay native
BNET models, rerun the full fixed-point census, or reconstruct full-model support
from native model files. Small displayed residual/component maps are evaluated
exhaustively. Retained full-model computation records are used as numerical inputs.
Run from any directory. Results default to ../qa/regulatory_verification.json.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re

FROZEN_SHA256 = {'data/bernoulli_failure_g_A.csv': '46b3f549c6777803421b20fd8b6fb1f991873614898bd09774f466e1cc08a200',
 'data/bernoulli_failure_g_B.csv': 'd19f4fba2199959cad0900681614bffcc5c94cd27b809bd0a8e87ef06af61de7',
 'data/bernoulli_failure_g_C.csv': '85d91917f1e5d69d46806a50b79a7872eabbe56472c868a52c0bffdbf6843b21',
 'data/bernoulli_failure_g_D.csv': '7fe7c93858d1be174108d75eea1dc61fdc38e5bbbea0aaedf3c799770ec0cf14',
 'data/bernoulli_failure_g_E.csv': '71e2f4348efdbf668be8b151389aa1f6640b433f0e8257ea40e7117d4cfdbf3f',
 'data/bernoulli_failure_g_prime_A.csv': '554b37e3e42e1794838ea53e5b71b87547a5b7ed3259bde2175432e553449fc6',
 'data/bernoulli_failure_g_prime_B.csv': '66809e4213def9e52932a0c94b22c6e7782abcb47c84431e55cc922efdb4d9c3',
 'data/bernoulli_failure_g_prime_C.csv': '0fe0e25dd47c6190463fe1e5082552bc5fe17f8fd3757c928172a46ca143a0c1',
 'data/bernoulli_failure_g_prime_D.csv': '46254d77ee43fb3a4e7b8ac32aefa88e8789941ede9da18598b40adf99bef1ae',
 'data/bernoulli_failure_g_prime_E.csv': '9b24c250301c75bad3b9ad464df53127c60f486f168a31b8641f1ec3177f7f23',
 'data/bernoulli_moments_A.csv': '7328c8cb446f30951fb62f51540cbd067f3d54a594059d4fdf75681d9d09a59e',
 'data/bernoulli_moments_B.csv': 'df15f545b742481a36effc7c4ac728f723f82a1a8954d28e4e5204997aca0c4f',
 'data/bernoulli_moments_C.csv': 'b942f830511862e6419ce165ed272b2dcd67bb80d01861f3b024076eea5a4dd0',
 'data/bernoulli_moments_D.csv': '36b06f0f755a6b15346078ebfdb4f54d60030c7637b195c52268aaf0c583cf71',
 'data/bernoulli_moments_E.csv': '73c501c041bace36b6b78c868a8427de5fb98ba36a22c9fd8db7a7f4a1cbb8e2',
 'data/bernoulli_pmf_A.csv': '5e69e28508407c5c78f472ecb1ccf8dcd609cd14e4b9c11b642f2576a6ae3224',
 'data/bernoulli_pmf_B.csv': '0cffcb05c1c5b25742f3d75036bbc083df2a7473311fb0b213917d2332b706e0',
 'data/bernoulli_pmf_C.csv': '81409d766efa92e4dd3990fe2a37368f8a94d448514271f1a338f8792db342a2',
 'data/bernoulli_pmf_D.csv': '38df854dd23255f8a1c4dda71e8937ab1a481cad707b85fb13066682c5bac795',
 'data/bernoulli_pmf_E.csv': 'be178a0162621755f5d52b68e42690a429da528f4e9af2a5e1eb8da3940523f1',
 'data/bernoulli_survival_A.csv': 'd9fcb246317701c1fb684c9b8697b39856b01af287dcb9fdc457961894924531',
 'data/bernoulli_survival_B.csv': 'ea34bc096abf6e70dd2dcfa960a06d06e07eb3f3b7be2f202ad64a245e4b1ed5',
 'data/bernoulli_survival_C.csv': '6cb16659e53d8d89593afc3beafed3b1f7bddb6ffd42c2547986040e2c8c9669',
 'data/bernoulli_survival_D.csv': '7adf824198b878ae2b1b691f3f246daaf1ff3cc900ba0296f4f817cf39b0d69d',
 'data/bernoulli_survival_E.csv': '247c115081cd3baab030f2f42877cccad7cbe8ede8062e121ebf7fea3a034e0a',
 'data/bernoulli_zero_g_A.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'data/bernoulli_zero_g_B.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'data/bernoulli_zero_g_C.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'data/bernoulli_zero_g_D.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'data/bernoulli_zero_g_E.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'data/bernoulli_zero_g_prime_A.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'data/bernoulli_zero_g_prime_B.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'data/bernoulli_zero_g_prime_C.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'data/bernoulli_zero_g_prime_D.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'data/bernoulli_zero_g_prime_E.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'data/fig4_convergence.csv': '73a3560b0ab23fdccfcf152a3a307259e32130e0784d1f6c757408c2a7c3e258',
 'data/fig5_cd4_apc.csv': 'cb2a376cfed7d2ebd6200a7ebe5f41b8404ffe8e4cfea8d667ec9bfc12d489fb',
 'data/fig5_mac_bacteria_constant_coordinates.csv': '03d74f07d872d7013571f7501455db9057c33686e90053076fbad1edc1dd7bb1',
 'data/fig5_mac_bacteria_cycle.csv': '885a27757741ca4e5fed84b4b6230e1662ca33cf4f98f887d4e897e83f9fb656',
 'data/fig5_mac_bacteria_full_orbit.csv': 'cb0b060e4e2690c3c0c1c157d2d4822fe9cb2c2e98bec917b247a44d0b1073ee',
 'data/fig5_mac_bacteria_status.csv': 'd4612a77c4c19a11646b1ad738e5e2e0ad518df6aa8dc77cbd3d4590527af2b3',
 'data/fig6_apc_isolated_rules.csv': '3bd9b43e50e3718e8363c4e9ea278b4d3e5516ec26b5185c01506e051ba0cb06',
 'data/fig6_bar_labels.tex': 'aa478e27d0108e7b2bb7d3bcfa01694fbaecb70150b29399d76122e8e260ae69',
 'data/fig6_build_metadata.json': '1aeb4cbf19f1ffdf4835d5467e1bc90d5b3a5fbda6f1a2a02809d452baefda66',
 'data/fig6_cd4_cost_distribution.csv': 'c285d1e2df095f640c3221f9a9d6052c3ab54d156c72d65a0ecbc955d9d1d596',
 'data/fig6_source_records.json': '1e89a9d147a00c56acfc5e6199214befcce28cbe5a67244473d8f1156d4b0863',
 'fig1_cost_differences.pdf': '45789c85b580a151b63618caddc03f11c1a8130d3bad17cc530bdc3f8290cb95',
 'fig1_cost_differences.tex': '9bf1b8142dd2ef122a3dd41c4bf064f79863d394c50eb517b3ce4cda41a2e5d9',
 'fig1_cost_distribution.pdf': '99a143b86b65edf847c9168d11e552e984fe77af1542b7d68a102a8af759549f',
 'fig1_cost_distribution.tex': '1ec68082b3759585050a9a1b9c8c86866d7a1077ed400207b50e9661681e44a4',
 'fig2_certificate_failure.pdf': '41b959362eef2fbcefea4cc70579fe352c95f01745c4d377f9829dda1f84e2b7',
 'fig2_certificate_failure.tex': '3b7a27d7dd5bc1a5c8b17518dc26e2b474151edc9ee5c930d4b0ea25b0a00f6f',
 'fig3_residual_and_moments.pdf': '4aa6fde81cccfec871f78ffbe48cc67412e1f84c953d456941e749d14785c186',
 'fig3_residual_and_moments.tex': 'ca6b9c372ddaae66f0ab701b291c56c1c7054f7a79ad34e4ab06aaec099eaa90',
 'fig_condition_moment_overview.pdf': '29ecc846417a9236eb838d168251c34c5b7c5404e233b421ed2b2d7d3f969579',
 'fig_condition_moment_overview.tex': '442d27953cbf3a27044ddb631a06bdc18a2dbc19596718feaa60b48d987a96d2'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def tex_without_comments(path):
    return re.sub(r'(?<!\\)%[^\n]*', '', path.read_text())


def table_mappings(source, filename):
    found = []
    for options, name in re.findall(r'table\s*\[([^]]*)\]\s*\{([^}]+)\}', source):
        if name == 'data/' + filename:
            x = re.search(r'(?:^|,)\s*x\s*=\s*([^,]+)', options)
            y = re.search(r'(?:^|,)\s*y\s*=\s*([^,]+)', options)
            require(x is not None and y is not None, 'Missing x/y table mapping')
            found.append((x.group(1).strip(), y.group(1).strip()))
    return found


def check_integrity(figures, baseline):
    for relative, expected in FROZEN_SHA256.items():
        require(sha(figures / relative) == expected, 'Retained asset changed: ' + relative)
        if baseline is not None:
            require(sha(figures / relative) == sha(baseline / relative),
                    'Baseline differs: ' + relative)
    expected_csv = {n for n in FROZEN_SHA256 if n.endswith('.csv')}
    actual_csv = {p.relative_to(figures).as_posix() for p in (figures / 'data').glob('*.csv')}
    require(actual_csv == expected_csv, 'Retained CSV inventory differs')
    return {'csv_files': len(expected_csv), 'retained_json_and_labels': 3,
            'unchanged_numerical_tex_and_pdf': 10, 'frozen_assets': len(FROZEN_SHA256),
            'optional_baseline_comparison': baseline is not None}


def check_cd4_support(figures):
    data = rows(figures / 'data/fig5_cd4_apc.csv')
    expected = [
        (0,154,150,0,0),(1,132,119,0,0),(2,117,102,0,0),(3,94,80,0,0),
        (4,75,60,0,0),(5,63,47,0,0),(6,49,32,0,0),(7,36,19,0,0),
        (8,26,9,0,0),(9,22,5,0,0),(10,20,3,0,0),(11,19,2,0,0),
        (12,18,1,0,0),(13,17,0,0,0),(14,17,0,1,1)]
    columns = ('support_round','baseline_support','optimal_support',
               'baseline_stationary_extension','optimal_stationary_extension')
    require([tuple(int(r[c]) for c in columns) for r in data] == expected,
            'CD4 support values or stationary flags differ')
    for name in ('fig5_larger_network_support.tex',):
        source = tex_without_comments(figures / name)
        require(table_mappings(source, 'fig5_cd4_apc.csv') ==
                [('support_round','baseline_support'),('support_round','optimal_support')],
                name + ': support columns differ')
    return {'rows': 15, 'initial_support': [154,150], 'terminal_support': [17,0],
            'first_empty_controlled_round': 13, 'stationary_extension_round': 14,
            'scope': 'Retained coordinate-cube support bounds; not exact reachable-image supports or a trajectory'}


def check_macrophage(figures):
    data = figures / 'data'
    orbit = rows(data / 'fig5_mac_bacteria_full_orbit.csv')
    names = list(orbit[0])[1:]
    require(len(names) == len(set(names)) == 321, 'Macrophage coordinate inventory differs')
    require([int(r['synchronous_update']) for r in orbit] == list(range(16)),
            'Macrophage retained update range differs')
    states = [tuple(int(r[n]) for n in names) for r in orbit]
    require(all(v in (0,1) for s in states for v in s), 'Nonbinary orbit entry')
    repetitions = [(i,j) for j in range(1,16) for i in range(j) if states[i] == states[j]]
    require(repetitions == [(11,15)], 'First full-state recurrence is not exactly (11,15)')
    varying = [i for i in range(321) if len({states[t][i] for t in range(11,15)}) > 1]
    cycle = rows(data / 'fig5_mac_bacteria_cycle.csv')
    constant = rows(data / 'fig5_mac_bacteria_constant_coordinates.csv')
    require(len(cycle) == 19 and len(constant) == 302, 'Cycle/constant partition size differs')
    require([int(r['node_index_zero_based']) for r in cycle] == varying,
            'Raster does not contain every varying coordinate in source order')
    expected_cells = {}
    for j, r in enumerate(cycle):
        i = int(r['node_index_zero_based'])
        require(r['source_node_name'] == names[i], 'Cycle coordinate identity mismatch')
        require(r['display_label'] == names[i].removeprefix('v_'), 'Unexpected display-name change')
        for t in range(11,16):
            require(int(r['update_' + str(t)]) == states[t][i], 'Orbit projection differs')
            expected_cells[(float(t),float(j))] = states[t][i]
    constant_indices = []
    inputs = []
    for r in constant:
        i = int(r['node_index_zero_based'])
        constant_indices.append(i)
        require(r['source_node_name'] == names[i], 'Constant coordinate identity mismatch')
        require(all(states[t][i] == int(r['value_on_cycle']) for t in range(11,16)),
                'Constant coordinate varies on retained cycle')
        require(int(r['external_input']) in (0,1), 'Invalid external-input flag')
        if int(r['external_input']):
            inputs.append(i)
    require(sorted(varying + constant_indices) == list(range(321)), 'Coordinate partition is not exact')
    require(len(inputs) == 19, 'Input count differs')
    one_inputs = {names[i] for i in inputs if states[0][i]}
    require(one_inputs == {'v_Bacteria','v_External_Activator'}, 'Clamped input context differs')
    require(all(states[t][i] == states[0][i] for i in inputs for t in range(16)),
            'External inputs are not clamped throughout retained orbit')
    require(all(states[0][i] == 0 for i in range(321) if i not in inputs),
            'Initial regulated state is not zero')
    status = rows(data / 'fig5_mac_bacteria_status.csv')
    require(len(status)==1 and status[0]['cycle_start_update']=='11'
            and status[0]['cycle_length']=='4'
            and status[0]['fixed_target_certificate_applied']=='False',
            'Retained orbit status differs')
    raster_checks = {}
    number = r'(-?\d+(?:\.\d+)?)'
    cell_pattern = (r'\\node\s*(?:\[[^]]*\])?\s*at\s*\(axis cs:\s*' + number +
                    r'\s*,\s*' + number + r'\s*\)\s*\{([01])\}\s*;')
    rectangle_pattern = (r'\\path\s*\[([^]]*)\]\s*\(axis cs:\s*' + number + r'\s*,\s*' + number +
                         r'\s*\)\s*rectangle\s*\(axis cs:\s*' + number + r'\s*,\s*' + number + r'\s*\)\s*;')
    for filename in ('fig5_larger_network_support.tex',):
        source = tex_without_comments(figures / filename)
        entries = re.findall(cell_pattern, source)
        rendered = {(float(x),float(y)):int(v) for x,y,v in entries}
        require(len(entries) == 95 and rendered == expected_cells, filename + ': 95 raster labels differ')
        rectangles = re.findall(rectangle_pattern, source)
        fills = {0:set(),1:set()}
        locations = []
        for options,x0,y0,x1,y1 in rectangles:
            x0,y0,x1,y1 = map(float,(x0,y0,x1,y1))
            loc = ((x0+x1)/2,(y0+y1)/2)
            require(x1-x0 == 1 and y1-y0 == 1 and loc in expected_cells,
                    filename + ': cell rectangle does not align with data bin')
            fill = re.search(r'(?:^|,)\s*fill\s*=\s*([^,]+)', options)
            require(fill is not None, filename + ': cell has no explicit fill')
            fills[expected_cells[loc]].add(fill.group(1).strip())
            locations.append(loc)
        require(len(locations)==95 and len(set(locations))==95, filename + ': raster rectangles differ')
        require(len(fills[0])==len(fills[1])==1 and fills[0]!=fills[1],
                filename + ': binary fills are not distinct and consistent')
        labels = re.search(r'yticklabels\s*=\s*\{(.*?)\},\s*\n', source, re.S)
        require(labels is not None, filename + ': no row labels')
        actual = re.findall(r'\{([^{}]*)\}', labels.group(1))
        require([s.replace(r'\_', '_') for s in actual] == [r['display_label'] for r in cycle],
                filename + ': raster row order/labels differ')
        raster_checks[filename] = {'cells':len(entries),'rows':19,'columns':5}
    return {'coordinates':321,'retained_states':16,'retained_adjacent_pairs':15,
            'first_repeat':[11,15],'period':4,'varying_cycle_coordinates':19,
            'constant_cycle_coordinates':302,'external_inputs':19,
            'zero_initial_regulated_coordinates':302,'raster_checks':raster_checks,
            'scope':'Checks retained full-state records and plotted projection; does not replay native BNET transitions or repeat exhaustive no-fixed-point enumeration'}


def minimum_pins(names, update):
    """Exhaust all permanent zero pins and every compatible initial state."""
    n = len(names)
    require(update(0)==0, 'The displayed residual target is not a fixed point')
    successful = []
    worst = {}
    initial_state_cases = 0
    for mask in range(1 << n):
        compatible = [s for s in range(1 << n) if not s & mask]
        good = True
        longest = 0
        for start in compatible:
            initial_state_cases += 1
            seen = set()
            state = start
            while state and state not in seen:
                seen.add(state)
                state = update(state) & ~mask
            if state:
                good = False
            else:
                longest = max(longest, len(seen))
        if good:
            successful.append(mask)
            worst[mask] = longest
    optimum = min(mask.bit_count() for mask in successful)
    minima = [mask for mask in successful if mask.bit_count()==optimum]
    return {'minimum':optimum,
            'minimum_sets':[[names[j] for j in range(n) if mask & (1 << j)] for mask in minima],
            'minimum_set_max_transients':[worst[m] for m in minima],
            'all_pin_sets_checked':1<<n,
            'compatible_initial_state_cases':initial_state_cases,
            'successful_pin_masks':successful}


def check_cd4_components(figures):
    data = figures / 'data'
    record = json.loads((data / 'fig6_source_records.json').read_text())
    require(record['schema'] == 'Regulatory-Component-Figure-v1',
            'Unexpected regulatory component record schema')
    require(record['source_json'] == 'workflows/biological/results/all_targets.json',
            'Unexpected retained regulatory source identity')
    distribution = rows(data / 'fig6_cd4_cost_distribution.csv')
    require(len(distribution)==6, 'CD4 context count differs')
    require([r['context'] for r in distribution] == record['context_order'], 'CD4 context order differs')
    total = 0
    target_ids = []
    expected_label_nodes = []
    for row, saved in zip(distribution, record['panel_a']):
        require(row['context']==saved['context'], 'Context identity differs')
        for key in ('row','cost2','cost3','cost4','cost5','total'):
            require(int(row[key])==saved[key], 'CD4 retained histogram differs at '+key)
        counts = Counter(t['exact_cost'] for t in saved['target_costs'])
        require(all(counts[k]==int(row['cost'+str(k)]) for k in range(2,6)), 'Target-by-target cost tally differs')
        require(sum(counts.values())==int(row['total']), 'Target total differs')
        target_ids.extend(t['target_id'] for t in saved['target_costs'])
        total += int(row['total'])
        start = 0
        for k in range(2,6):
            count = int(row['cost'+str(k)])
            if count:
                expected_label_nodes.append((start+count/2,float(row['row']),count))
            start += count
        expected_label_nodes.append((float(start),float(row['row']),start))
    require(total==46 and len(set(target_ids))==46, 'CD4 census does not contain 46 unique records')
    labels = (data / 'fig6_bar_labels.tex').read_text()
    actual_labels = [(float(x),float(y),int(v)) for x,y,v in re.findall(
        r'\(axis cs:([\d.]+),([\d.]+)\)\s*\{(\d+)\}', labels)]
    require(actual_labels==expected_label_nodes, 'Stacked-bar segment or total label geometry differs')
    source = tex_without_comments(figures / 'fig6_cd4_target_mechanisms.tex')
    require(table_mappings(source, 'fig6_cd4_cost_distribution.csv')==
            [('cost'+str(k),'row') for k in range(2,6)], 'CD4 stacked-bar mappings differ')
    require('xbar stacked' in source and r'\input{data/fig6_bar_labels.tex}' in source,
            'CD4 stacked-bar geometry or retained labels missing')
    normalized = re.sub(r'\s+', '', source)
    normalized = re.sub(r'\\mathrm\{([^{}]+)\}', r'\1', normalized)
    normalized = normalized.replace('^{+}', '^+').replace(r'\qquad','').replace(r'\quad','')
    require(Counter(re.findall(r'\\kappa_B\(C\)=(\d+)', normalized)) == Counter({'1':4,'0':1}),
            'Displayed isolated-component costs differ')
    require('c_B=1+1+1+1+0=4' in normalized, 'Displayed APC cost sum differs')
    expected_expressions = {'Dec2':'GATA3','GATA3':'Dec2','FOXP3':r'FOXP3\veeSMAD3',
                            'SMAD3':'TGFBR','TGFB':'0','TGFBR':'TGFB','RORGT':'0',
                            'TBET':'TBET','p':'p'}
    symbolic = {
        'v_Dec2':lambda s:s['v_GATA3'], 'v_GATA3':lambda s:s['v_Dec2'],
        'v_FOXP3':lambda s:s['v_FOXP3'] | s['v_SMAD3'],
        'v_SMAD3':lambda s:s['v_TGFBR'], 'v_TGFB':lambda s:0,
        'v_TGFBR':lambda s:s['v_TGFB'], 'v_RORGT':lambda s:0,
        'v_TBET':lambda s:s['v_TBET'], 'v_proliferation':lambda s:s['v_proliferation']}
    for name, expression in expected_expressions.items():
        matches = re.findall(re.escape(name+'^+=')+r'(.*?)(?=\$|,|;)', normalized)
        require(matches == [expression], 'Displayed isolated rule differs: '+name)
    isolated = rows(data / 'fig6_apc_isolated_rules.csv')
    require(len(isolated)==9, 'Isolated rule count differs')
    panel = record['panel_b']
    require(len(panel['target_by_name'])==188, 'Selected CD4 target size differs')
    require(all(panel['target_by_name'][n]==0 for n in symbolic), 'Displayed target coordinates are not zero')
    require(panel['target_id']=='CD4_APC_FP01', 'Selected target identity differs')
    require(len(panel['acyclic_singleton_names'])==8 and len(panel['residual_names'])==17,
            'Residual component partition size differs')
    require(set(panel['residual_names']) == set(symbolic) | set(panel['acyclic_singleton_names']),
            'Residual component partition differs')
    truth_values = 0
    components = []
    expected_minima = [[['v_Dec2'],['v_GATA3']],[['v_FOXP3']],[[]],[['v_TBET']],[['v_proliferation']]]
    for j, component in enumerate(panel['cyclic_components']):
        names = component['names']
        tables = {}
        csv_rows = [r for r in isolated if int(r['component'])==j]
        require([r['node'] for r in csv_rows]==names, 'Component row order differs')
        for saved, row in zip(component['isolated_target_error_rows'],csv_rows):
            parents = row['parents_little_endian'].split(';') if row['parents_little_endian'] else []
            table = [int(v) for v in row['table']]
            require(saved=={'name':row['node'],'parents':parents,'table':table}, 'CSV/source-record isolated rule mismatch')
            require(int(row['component_cost'])==component['cost'], 'Component cost metadata differs')
            require(len(table)==1<<len(parents), 'Incomplete local truth table')
            tables[row['node']] = (parents,table)
        for state in range(1<<len(names)):
            assignment = {name:(state>>i)&1 for i,name in enumerate(names)}
            for name in names:
                parents, table = tables[name]
                index = sum(assignment[p]<<i for i,p in enumerate(parents))
                require(table[index]==symbolic[name](assignment), 'Symbolic/truth-table output mismatch')
                truth_values += 1
        def update(state):
            assignment = {name:(state>>i)&1 for i,name in enumerate(names)}
            return sum(symbolic[name](assignment)<<i for i,name in enumerate(names))
        result = minimum_pins(names, update)
        require(result['minimum']==component['cost'] and result['minimum_sets']==expected_minima[j],
                'Exact isolated component optimum differs')
        result.pop('successful_pin_masks')
        components.append({'names':names,**result})
    require(truth_values==78 and panel['independently_checked_source_rule_values']==78,
            'Complete-component output count differs')
    require(sum(c['minimum'] for c in components)==panel['exact_additional_cost']==4,
            'Selected APC cost decomposition differs')
    return {'target_records':total,'contexts':6,'bar_labels':len(actual_labels),
            'isolated_rules':9,'full_component_truth_values':truth_values,'components':components,
            'acyclic_zero_cost_singletons':8,'selected_apc_exact_additional_cost':4,
            'scope':'Retained exact target-cost tally; exhaustive evaluation of displayed isolated maps. Does not replay the 188-coordinate native model or recompute full SCCs/census.'}


def check_small_certificates(figures):
    source = tex_without_comments(figures / 'fig4_biological_certificates.tex')
    arrows = []
    for options, body in re.findall(r'\\draw\s*\[([^]]*)\]\s*(.*?);',source,re.S):
        vertices = re.findall(r'\(v(\d+)\)',body)
        if vertices and ('->' in options or '<-' in options):
            require(len(vertices)==2 and '->' in options and '<-' not in options,
                    'Copy schematic has an unexpected undirected/reversed edge')
            arrows.append(tuple(map(int,vertices)))
    expected = [(15,16),(16,15),(15,17),(17,18),(18,19)]
    require(Counter(arrows)==Counter(expected),'T-LGL copy-arrow topology differs')
    for node in range(15,20):
        require(re.search(r'\(v'+str(node)+r'\).*?\{\$y_\{'+str(node)+r'\}\$\}',source,re.S),
                'Copy schematic node label mismatch')
    parents = [1,0,0,2,3]
    copy_update = lambda state:sum(((state>>p)&1)<<j for j,p in enumerate(parents))
    copy_result = minimum_pins(list(range(15,20)),copy_update)
    require(copy_result['minimum']==1 and copy_result['minimum_sets']==[[15],[16]],
            'Copy graph minimum-pin condition differs')
    require(all(bool(mask&3)==(mask in copy_result['successful_pin_masks']) for mask in range(32)),
            'Copy graph convergence does not equal hitting its cycle')
    copy_result.pop('successful_pin_masks')
    data = rows(figures / 'data/fig4_convergence.csv')
    expected_rows = [(0,14,14,16384),(1,12,12,720),(2,11,11,135),(3,9,9,56),
                     (4,7,7,24),(5,7,7,12),(6,7,7,7),(7,7,6,5),(8,7,5,4),
                     (9,7,3,3),(10,7,1,2),(11,7,0,1)]
    cols = ('update','cube_support_size','exact_reachable_union_support_size','distinct_reachable_states')
    require([tuple(int(r[k]) for k in cols) for r in data]==expected_rows,
            'Segment-polarity retained supports/reachable counts differ')
    require(table_mappings(source,'fig4_convergence.csv')==
            [('update','cube_support_size'),('update','exact_reachable_union_support_size')],
            'Segment-polarity support mapping differs')
    def residual_update(state):
        a,b,c,d,e,f,h = [(state>>i)&1 for i in range(7)]
        values = (h & (a | f), a, f & h, d & (1-c), d, e, e)
        return sum(v<<i for i,v in enumerate(values))
    times = []
    for initial in range(128):
        state,seen = initial,set()
        while state and state not in seen:
            seen.add(state)
            state = residual_update(state)
        require(state==0,'Displayed seven-coordinate residual does not converge')
        times.append(len(seen))
    require(max(times)==8 and times[8]==8, 'Residual maximum transient/witness differs')
    content_path = figures / 'FIGURE_CONTENT.json'
    content = json.loads(content_path.read_text())
    expected_residual = {'coordinate_order': ['a', 'b', 'c', 'd', 'e', 'f', 'h'],
     'node_indices_one_based': [4, 5, 10, 11, 13, 16, 17],
     'update_expressions': {'a': 'h & (a | f)',
                            'b': 'a',
                            'c': 'f & h',
                            'd': 'd & (1-c)',
                            'e': 'd',
                            'f': 'e',
                            'h': 'e'},
     'target': [0, 0, 0, 0, 0, 0, 0],
     'update_convention': 'synchronous'}
    require(content['segment_residual'] == expected_residual,
            'Retained segment-polarity residual map differs from the evaluated map')
    return {'copy_arrows':5,'copy_graph':copy_result,
            'segment_retained_rows':12,'full_compatible_initial_states':16384,
            'full_maximum_transient_from_retained_support':11,'stationary_cube_support':7,
            'residual_states_exhausted':128,'residual_maximum_transient':8,
            'retained_residual_map_matched':True,
            'retained_residual_map_source':'FIGURE_CONTENT.json/segment_residual',
            'retained_residual_map_record_sha256':sha(content_path),
            'residual_eight_update_witness':[0,0,0,1,0,0,0],
            'scope':'Full-system 11-update claim checked against retained exhaustive-support data; seven-coordinate residual independently enumerated. Native full-model replay is not performed.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figures',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--baseline-figures',type=Path,
                        help='Optional original figures directory for extra byte-preservation comparison')
    parser.add_argument('--output',type=Path,help='Result JSON (default: source-root qa/regulatory_verification.json)')
    args = parser.parse_args()
    figures = args.figures.resolve()
    report = {'status':'PASS','scope':__doc__.split('\n\n')[1].strip(),
              'dependency':'Python standard library only','checks':{}}
    checks = [('integrity',lambda:check_integrity(figures,args.baseline_figures)),
              ('cd4_support',lambda:check_cd4_support(figures)),
              ('macrophage_orbit_and_raster',lambda:check_macrophage(figures)),
              ('cd4_target_costs_and_components',lambda:check_cd4_components(figures)),
              ('small_residual_certificates',lambda:check_small_certificates(figures))]
    for label,check in checks:
        try:
            report['checks'][label] = {'status':'PASS',**check()}
        except (ValueError,OSError,KeyError,IndexError,TypeError) as exc:
            report['status'] = 'FAIL'
            report['checks'][label] = {'status':'FAIL','reason':str(exc)}
    output = args.output or figures.parent / 'qa/regulatory_verification.json'
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':report['status'],'checks':{k:v['status'] for k,v in report['checks'].items()},
                      'result':str(output)},indent=2))
    raise SystemExit(0 if report['status']=='PASS' else 1)


if __name__ == '__main__':
    main()
