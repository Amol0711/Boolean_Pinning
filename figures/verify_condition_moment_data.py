#!/usr/bin/env python3
"""Check retained condition-failure, residual-survival and moment figures.

Python 3.10 or newer is required. SciPy is optional: when available, the script
recalculates all exact binomial endpoints; otherwise it explicitly reports that
recalculation as skipped. Included aggregate hashes, arithmetic and source
mappings are checked in either case. This is not a raw-network replay and does
not recompute the approximate moment standard errors.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re

try:
    from scipy.stats import beta
except ImportError:
    beta = None

M = 10000
SIZES = [100, 400, 1600, 6400, 25600]
OFFSETS = {'A':.86, 'B':.93, 'C':1., 'D':1/.93, 'E':1/.86}
TYPES = {'A':((1,.5),(3,.1)), 'B':((1,.1),(3,.5)),
         'C':((1,.8),(3,.2)), 'D':((1,.35),(3,.35)), 'E':((0,.8),(4,.3))}
FAILURE_COUNTS = {
    'g': {'A':[178,51,8,1,1], 'B':[4219,3053,1839,810,274],
          'C':[1716,771,294,71,20], 'D':[2681,1481,522,182,52], 'E':[2642,1417,533,166,46]},
    'g_prime': {'A':[115,31,6,1,1], 'B':[3024,1641,664,201,48],
                'C':[990,341,111,32,6], 'D':[1773,707,187,39,13], 'E':[1531,568,152,39,9]}}
MAXIMA = {'A':15,'B':199,'C':105,'D':109,'E':125}
LARGEST_MOMENTS = {'A':[[.5149,.5149],[.7823,.7823]],
                   'B':[[1.5863,1.5870],[4.1035,4.1070]],
                   'C':[[1.2031,1.2031],[2.6461,2.6461]],
                   'D':[[1.1890,1.1890],[2.6026,2.6026]],
                   'E':[[.9252,.9252],[1.7722,1.7722]]}

CSV_SHA256 = {'bernoulli_failure_g_A.csv': '46b3f549c6777803421b20fd8b6fb1f991873614898bd09774f466e1cc08a200',
 'bernoulli_failure_g_B.csv': 'd19f4fba2199959cad0900681614bffcc5c94cd27b809bd0a8e87ef06af61de7',
 'bernoulli_failure_g_C.csv': '85d91917f1e5d69d46806a50b79a7872eabbe56472c868a52c0bffdbf6843b21',
 'bernoulli_failure_g_D.csv': '7fe7c93858d1be174108d75eea1dc61fdc38e5bbbea0aaedf3c799770ec0cf14',
 'bernoulli_failure_g_E.csv': '71e2f4348efdbf668be8b151389aa1f6640b433f0e8257ea40e7117d4cfdbf3f',
 'bernoulli_failure_g_prime_A.csv': '554b37e3e42e1794838ea53e5b71b87547a5b7ed3259bde2175432e553449fc6',
 'bernoulli_failure_g_prime_B.csv': '66809e4213def9e52932a0c94b22c6e7782abcb47c84431e55cc922efdb4d9c3',
 'bernoulli_failure_g_prime_C.csv': '0fe0e25dd47c6190463fe1e5082552bc5fe17f8fd3757c928172a46ca143a0c1',
 'bernoulli_failure_g_prime_D.csv': '46254d77ee43fb3a4e7b8ac32aefa88e8789941ede9da18598b40adf99bef1ae',
 'bernoulli_failure_g_prime_E.csv': '9b24c250301c75bad3b9ad464df53127c60f486f168a31b8641f1ec3177f7f23',
 'bernoulli_moments_A.csv': '7328c8cb446f30951fb62f51540cbd067f3d54a594059d4fdf75681d9d09a59e',
 'bernoulli_moments_B.csv': 'df15f545b742481a36effc7c4ac728f723f82a1a8954d28e4e5204997aca0c4f',
 'bernoulli_moments_C.csv': 'b942f830511862e6419ce165ed272b2dcd67bb80d01861f3b024076eea5a4dd0',
 'bernoulli_moments_D.csv': '36b06f0f755a6b15346078ebfdb4f54d60030c7637b195c52268aaf0c583cf71',
 'bernoulli_moments_E.csv': '73c501c041bace36b6b78c868a8427de5fb98ba36a22c9fd8db7a7f4a1cbb8e2',
 'bernoulli_pmf_A.csv': '5e69e28508407c5c78f472ecb1ccf8dcd609cd14e4b9c11b642f2576a6ae3224',
 'bernoulli_pmf_B.csv': '0cffcb05c1c5b25742f3d75036bbc083df2a7473311fb0b213917d2332b706e0',
 'bernoulli_pmf_C.csv': '81409d766efa92e4dd3990fe2a37368f8a94d448514271f1a338f8792db342a2',
 'bernoulli_pmf_D.csv': '38df854dd23255f8a1c4dda71e8937ab1a481cad707b85fb13066682c5bac795',
 'bernoulli_pmf_E.csv': 'be178a0162621755f5d52b68e42690a429da528f4e9af2a5e1eb8da3940523f1',
 'bernoulli_survival_A.csv': 'd9fcb246317701c1fb684c9b8697b39856b01af287dcb9fdc457961894924531',
 'bernoulli_survival_B.csv': 'ea34bc096abf6e70dd2dcfa960a06d06e07eb3f3b7be2f202ad64a245e4b1ed5',
 'bernoulli_survival_C.csv': '6cb16659e53d8d89593afc3beafed3b1f7bddb6ffd42c2547986040e2c8c9669',
 'bernoulli_survival_D.csv': '7adf824198b878ae2b1b691f3f246daaf1ff3cc900ba0296f4f817cf39b0d69d',
 'bernoulli_survival_E.csv': '247c115081cd3baab030f2f42877cccad7cbe8ede8062e121ebf7fea3a034e0a',
 'bernoulli_zero_g_A.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'bernoulli_zero_g_B.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'bernoulli_zero_g_C.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'bernoulli_zero_g_D.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'bernoulli_zero_g_E.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'bernoulli_zero_g_prime_A.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'bernoulli_zero_g_prime_B.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'bernoulli_zero_g_prime_C.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'bernoulli_zero_g_prime_D.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'bernoulli_zero_g_prime_E.csv': '37a9fbce99adf6a8fc8159933fb534fb9edc9fbf7c02e9bdf24fbaf9476e5ffd',
 'fig4_convergence.csv': '73a3560b0ab23fdccfcf152a3a307259e32130e0784d1f6c757408c2a7c3e258',
 'fig5_cd4_apc.csv': 'cb2a376cfed7d2ebd6200a7ebe5f41b8404ffe8e4cfea8d667ec9bfc12d489fb',
 'fig5_mac_bacteria_constant_coordinates.csv': '03d74f07d872d7013571f7501455db9057c33686e90053076fbad1edc1dd7bb1',
 'fig5_mac_bacteria_cycle.csv': '885a27757741ca4e5fed84b4b6230e1662ca33cf4f98f887d4e897e83f9fb656',
 'fig5_mac_bacteria_full_orbit.csv': 'cb0b060e4e2690c3c0c1c157d2d4822fe9cb2c2e98bec917b247a44d0b1073ee',
 'fig5_mac_bacteria_status.csv': 'd4612a77c4c19a11646b1ad738e5e2e0ad518df6aa8dc77cbd3d4590527af2b3',
 'fig6_apc_isolated_rules.csv': '3bd9b43e50e3718e8363c4e9ea278b4d3e5516ec26b5185c01506e051ba0cb06',
 'fig6_cd4_cost_distribution.csv': 'c285d1e2df095f640c3221f9a9d6052c3ab54d156c72d65a0ecbc955d9d1d596'}

PROTECTED_FIGURE_SHA256 = {'fig1_cost_differences.pdf': '45789c85b580a151b63618caddc03f11c1a8130d3bad17cc530bdc3f8290cb95',
 'fig1_cost_differences.tex': '9bf1b8142dd2ef122a3dd41c4bf064f79863d394c50eb517b3ce4cda41a2e5d9',
 'fig1_cost_distribution.pdf': '99a143b86b65edf847c9168d11e552e984fe77af1542b7d68a102a8af759549f',
 'fig1_cost_distribution.tex': '1ec68082b3759585050a9a1b9c8c86866d7a1077ed400207b50e9661681e44a4'}

MARKS = {'A':'*','B':'square*','C':'triangle*','D':'diamond*','E':'pentagon*'}
COLORS = {'A':'0072B2','B':'D55E00','C':'009E73','D':'CC79A7','E':'806400'}

def require(condition, explanation):
    if not condition:
        raise AssertionError(explanation)

def close(a, b, label, tolerance=2e-12):
    require(math.isfinite(a) and math.isfinite(b), label+': nonfinite value')
    require(abs(a-b) <= tolerance, f'{label}: {a} != {b}')

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read_rows(path):
    with path.open(newline='') as stream:
        return [{key:float(value) for key,value in row.items()} for row in csv.DictReader(stream)]

def check_data(source):
    data_dir = source/'figures/data'
    metadata = json.loads((source/'figures/bernoulli_plot_metadata.json').read_text())
    require(metadata['network_size_display_offsets'] == OFFSETS, 'Display offsets changed')
    require(metadata['sample_accounting'] == {'cells':25,'draws_per_cell':M,'total_draws':250000,'excluded_draws':0},
            'Complete-sample accounting changed')
    require(set(CSV_SHA256) == {p.name for p in data_dir.glob('*.csv')}, 'CSV inventory changed')
    hashes = {}
    for name, expected in CSV_SHA256.items():
        hashes[name] = digest(data_dir/name)
        require(hashes[name] == expected, 'Retained CSV integrity failed: '+name)
    failure_results, moments_results, survival_results = {}, {}, {}
    failure_rows = moment_rows = moment_ranges = survival_rows = zero_survival = 0
    max_cp_error = max_integer_grid_error = 0.
    moment_extrema = {1:[math.inf,-math.inf],2:[math.inf,-math.inf]}
    failure_extrema = [math.inf,-math.inf]
    observed_zero_failures = []
    for condition in ('g','g_prime'):
        failure_results[condition] = {}
        for mixture in TYPES:
            rows = read_rows(data_dir/f'bernoulli_failure_{condition}_{mixture}.csv')
            zero_rows = read_rows(data_dir/f'bernoulli_zero_{condition}_{mixture}.csv')
            require([int(row['n']) for row in rows] == SIZES, 'Incomplete failure size grid')
            require([int(row['count']) for row in rows] == FAILURE_COUNTS[condition][mixture], 'Changed failure counts')
            expected_zeros = []
            for row in rows:
                count = int(row['count'])
                close(row['count'], count, 'Integer failure count')
                frequency = count/M
                cp_lo = row['cp_lo'] if beta is None else (0. if count == 0 else beta.ppf(.025,count,M-count+1))
                cp_hi = row['cp_hi'] if beta is None else (1. if count == M else beta.ppf(.975,count+1,M-count))
                require(0 <= cp_lo <= frequency <= cp_hi <= 1, 'Invalid retained CP envelope')
                close(row['true_rate'], frequency, 'Full-sample failure denominator')
                close(row['x'], row['n']*OFFSETS[mixture], 'Failure display offset', 1e-9)
                close(row['cp_lo'], cp_lo, 'Independent CP lower endpoint')
                close(row['cp_hi'], cp_hi, 'Independent CP upper endpoint')
                max_cp_error = max(max_cp_error,abs(row['cp_lo']-cp_lo),abs(row['cp_hi']-cp_hi))
                close(row['minus'], frequency-cp_lo, 'Failure lower error length')
                close(row['plus'], cp_hi-frequency, 'Failure upper error length')
                if count:
                    close(row['p'], frequency, 'Positive observed failure frequency')
                    failure_extrema[0] = min(failure_extrema[0],cp_lo)
                    failure_extrema[1] = max(failure_extrema[1],cp_hi)
                else:
                    require(math.isnan(row['p']) and row['cp_lo'] == 0, 'Zero frequency is plotted as positive')
                    expected_zeros.append({'n':row['n'],'x':row['x'],'upper':cp_hi})
                    observed_zero_failures.append([condition,mixture,int(row['n'])])
                failure_rows += 1
            require(len(zero_rows) == len(expected_zeros), 'Zero-upper-limit inventory disagrees')
            for row,expected in zip(zero_rows,expected_zeros):
                require(row.keys() == expected.keys(), 'Unexpected zero-limit fields')
                for key in row:
                    close(row[key],expected[key], 'Zero-limit '+key)
            failure_results[condition][mixture] = [int(row['count']) for row in rows]
    for mixture in TYPES:
        require(all(a <= b for a,b in zip(failure_results['g_prime'][mixture],failure_results['g'][mixture])),
                'Componentwise copy condition failures exceed residual copy failures')
        lam = sum(degree*q for degree,q in TYPES[mixture])/2
        mu = -math.log1p(-lam)
        normalizers = {1:mu,2:mu+mu*mu}
        rows = read_rows(data_dir/f'bernoulli_moments_{mixture}.csv')
        require([int(row['n']) for row in rows] == SIZES, 'Incomplete moment size grid')
        unnormalized = []
        for row in rows:
            close(row['x'],row['n']*OFFSETS[mixture],'Moment display offset',1e-9)
            raw_bounds = {}
            for power in (1,2):
                lo,hi,outerlo,outerhi = [row[stem+str(power)] for stem in ('lo','hi','outerlo','outerhi')]
                require(0 <= outerlo <= lo <= hi <= outerhi, 'Moment range nesting failed')
                close(row['err'+str(power)],hi-lo,'Moment computational range')
                close(row['outerminus'+str(power)],lo-outerlo,'Moment lower sampling length')
                close(row['outerplus'+str(power)],outerhi-lo,'Moment upper sampling length from lower marker')
                raw_bounds[power] = [lo*normalizers[power],hi*normalizers[power]]
                for value in raw_bounds[power]:
                    error = abs(M*value-round(M*value))
                    max_integer_grid_error = max(max_integer_grid_error,error)
                    require(error < 1e-7, 'Inverse-normalized cost powers are not complete-sample integer sums')
                moment_extrema[power][0] = min(moment_extrema[power][0],outerlo)
                moment_extrema[power][1] = max(moment_extrema[power][1],outerhi)
                moment_ranges += 1
            for endpoint in (0,1):
                require(raw_bounds[2][endpoint]+1e-10 >= raw_bounds[1][endpoint]**2,
                        'Endpoint first/second moment Jensen inequality fails')
            unnormalized.append({'n':int(row['n']),'first_raw':raw_bounds[1],'second_raw':raw_bounds[2]})
            moment_rows += 1
        for power in (1,2):
            for actual,expected in zip(unnormalized[-1]['first_raw' if power == 1 else 'second_raw'],
                                       LARGEST_MOMENTS[mixture][power-1]):
                close(actual,expected,'Largest-size raw-moment consistency')
        moments_results[mixture] = {'lambda_star':lam,'normalizers':normalizers,
                                    'inverse_normalized_computational_bounds':unnormalized}
        rows = read_rows(data_dir/f'bernoulli_survival_{mixture}.csv')
        require([row['x'] for row in rows] == list(range(MAXIMA[mixture]+2)), 'Incomplete retained survival thresholds')
        close(rows[0]['survival'],1.,'Survival at threshold zero')
        require(all(a['survival'] >= b['survival'] for a,b in zip(rows,rows[1:])), 'Survival is not nonincreasing')
        counts=[]
        for row in rows:
            survival = row['survival']
            require(0 <= survival <= 1,'Invalid survival probability')
            close(M*survival,round(M*survival),'Integer survival count',1e-8)
            counts.append(round(M*survival))
            if survival:
                close(row['plot_survival'],survival,'Positive retained survival ordinate')
            else:
                require(math.isnan(row['plot_survival']),'Zero survival must be omitted only from the log display')
                zero_survival += 1
            survival_rows += 1
        require(counts[-1] == 0 and counts[-2] > 0, 'Observed maximum does not match terminal records')
        masses=[a-b for a,b in zip(counts,counts[1:])]
        require(sum(masses) == M and all(x >= 0 for x in masses), 'Survival does not define a complete empirical residual histogram')
        survival_results[mixture] = {'rows':len(rows),'observed_maximum':MAXIMA[mixture],
            'terminal_zero_threshold':MAXIMA[mixture]+1,'survival_counts':counts,'residual_histogram_counts':masses}
    require((failure_rows,moment_rows,moment_ranges,survival_rows,zero_survival) == (50,25,50,563,5),
            'Unexpected aggregate record coverage')
    for power in (1,2):
        for actual,expected in zip(moment_extrema[power],metadata['moment_limits'][str(power)]):
            close(actual,expected,'Retained full moment-envelope extent')
    require(metadata['residual_maxima'] == MAXIMA,'Retained residual maxima disagree')
    require(not observed_zero_failures and not metadata['zero_failure_cells'],'Unexpected observed zero-failure cell')
    unchanged_figures = {}
    for name, expected in PROTECTED_FIGURE_SHA256.items():
        unchanged_figures[name] = digest(source/'figures'/name)
        require(unchanged_figures[name] == expected, 'Protected probability figure changed: '+name)
    return {'status':'PASS' if beta is not None else 'PASS_WITH_SKIPPED_CP','scope':'Retained aggregate arithmetic; no raw-network replay or recomputation of moment standard errors',
            'failure_rows_checked':failure_rows,'cp_endpoints_recomputed':2*failure_rows if beta is not None else 0,
            'cp_recalculation':'PASS' if beta is not None else 'SKIPPED: SciPy is unavailable; retained hashes and arithmetic were checked',
            'moment_rows_checked':moment_rows,'normalized_moment_ranges_checked':moment_ranges,
            'survival_rows_checked':survival_rows,'zero_survival_rows_retained':zero_survival,
            'failure_counts':failure_results,'failure_totals':{k:sum(sum(v) for v in groups.values()) for k,groups in failure_results.items()},
            'maximum_cp_endpoint_difference':max_cp_error if beta is not None else None,'failure_outer_endpoint_extrema':failure_extrema,
            'maximum_inverse_normalized_integer_grid_error':max_integer_grid_error,
            'normalized_moment_outer_endpoint_extrema':moment_extrema,'moments':moments_results,'survival':survival_results,
            'all_csvs_byte_identical':True,'csv_files_checked':len(hashes),'csv_sha256':hashes,
            'unchanged_probability_figure_sha256':unchanged_figures,
            'observed_zero_failure_cells':observed_zero_failures,
            'limitations':['Approximate moment sampling envelopes are preserved and checked algebraically; their endpoint standard errors require raw endpoint samples or equivalent sufficient statistics.',
                           'CP intervals are pointwise, not simultaneous, and do not bound Poisson approximation error.',
                           'This figure-export check does not consume independent-statistics or raw-network records retained elsewhere in the repository.']}


def options_values(text):
    values = {key:float(value) for key,value in re.findall(
        r'\b(xmin|xmax|ymin|ymax)\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)',text)}
    for key,value in re.findall(r'\b(xmode|ymode)\s*=\s*(log|linear)',text):
        values[key]=value
    return values


def plot_commands(text):
    pattern=r'\\addplot\+?\s*\[([^\]]*)\]\s*table\s*\[([^\]]*)\]\s*\{([^}]+)\}\s*;'
    commands=[]
    for options,table,name in re.findall(pattern,text,re.S):
        columns=dict(re.findall(r'(x|y|y error plus|y error minus|y error)\s*=\s*([^,\]]+)',table))
        columns={k:v.strip() for k,v in columns.items()}
        require('restrict' not in options+table and 'filter' not in options+table,'Data filtering is not permitted')
        commands.append((re.sub(r'\s+','',options),columns,name))
    return commands


def inspect_figure(path):
    text=path.read_text()
    actual_colors=dict(re.findall(r'\\definecolor\{mix([A-E])\}\{HTML\}\{([^}]+)\}',text))
    require(actual_colors == COLORS,'Changed mixture color identities')
    panels=[]
    group_sizes=[]
    for group in re.findall(r'\\begin\{groupplot\}(.*?)\\end\{groupplot\}',text,re.S):
        chunks=group.split('\\nextgroupplot')
        group_size=re.search(r'group size\s*=\s*(\d+) by (\d+)',chunks[0])
        require(group_size is not None,'Missing group layout')
        group_sizes.append([int(group_size.group(1)),int(group_size.group(2))])
        common=options_values(chunks[0])
        for chunk in chunks[1:]:
            axis_options=chunk.split('\\addplot',1)[0]
            bounds=dict(common,**options_values(axis_options))
            bounds.setdefault('xmode','linear');bounds.setdefault('ymode','linear')
            commands=plot_commands(chunk)
            require(commands,'A scientific panel has no retained-data plot')
            files={name for options,columns,name in commands}
            first=commands[0][2]
            if 'bernoulli_failure_g_prime_' in first: role='g_prime'
            elif 'bernoulli_failure_g_' in first: role='g'
            elif 'bernoulli_survival_' in first: role='survival'
            elif 'bernoulli_moments_' in first:
                require(commands[0][1].get('y') in ('lo1','lo2'),'Unrecognized moment ordinate')
                role='moment'+commands[0][1]['y'][-1]
            else: raise AssertionError('Unrecognized plotted data: '+first)
            expected_bounds={
                'g':(70,36000,.0000015,.7,'log','log'),
                'g_prime':(70,36000,.0000015,.7,'log','log'),
                'survival':(0,200,.00007,1.2,'linear','log'),
                'moment1':(70,37000,.8275195516154046,1.2107929637406762,'log','linear'),
                'moment2':(70,37000,.5315788186407935,3.289914264413321,'log','linear')}[role]
            for key,value in zip(('xmin','xmax','ymin','ymax','xmode','ymode'),expected_bounds):
                if isinstance(value,str): require(bounds[key]==value,'Changed axis scale: '+role)
                else: close(bounds[key],value,'Preserved full axis '+key)
            if role!='survival':
                ticks=re.findall(r'xtick=\{([^}]+)\}',chunks[0]+axis_options)
                require(ticks and [int(v) for v in ticks[-1].split(',')]==SIZES,
                        'All five true size locations must be explicit ticks')
            functions=re.findall(r'\\addplot\+?\s*\[[^\]]*\]\s*\{([^}]+)\}\s*;',chunk)
            require(functions==(['1/x'] if role in ('g','g_prime') else ['1'] if role.startswith('moment') else []),
                    'Unexpected reference function or reference attached to wrong panel')
            if role.startswith('moment'):
                power=role[-1]
                require('\\widehat m_{'+power+',n}' in axis_options,'Missing size-indexed moment notation')
                require('\\mu+\\mu^2' in axis_options if power=='2' else '/\\mu' in axis_options,
                        'Wrong raw-moment normalization label')
            elif role=='survival':
                require('|R|\\geq r' in axis_options,'Incorrect residual survival event')
            else:
                condition="G'_n" if role=='g_prime' else 'G_n'
                require(condition in axis_options,'Condition-specific title changed')
                require('Observed failure frequency' in text or '\\widehat{\\mathbb P}' in axis_options,
                        'Observed condition-failure meaning is not stated')
            observed={}
            endpoint_count=0
            for options,columns,name in commands:
                match=re.fullmatch(r'data/bernoulli_(failure_(g|g_prime)|survival|moments)_([A-E])\.csv',name)
                require(match is not None,'Unexpected source table: '+name)
                mixture=match.group(3)
                require('color=mix'+mixture in options,'Changed mixture color')
                require(columns.get('x')=='x','True/display abscissa confused')
                rows=read_rows(path.parent/name)
                if role in ('g','g_prime'):
                    require(match.group(2)==role and columns.get('y')=='p'
                            and columns.get('y error plus')=='plus' and columns.get('y error minus')=='minus',
                            'Wrong condition or probability error columns')
                    require('ydir=both' in options,'Failure interval must have both endpoints')
                    kind='failure'; endpoint_values=lambda row:[row['cp_lo'],row['cp_hi']]
                elif role=='survival':
                    require(match.group(1)=='survival' and columns.get('y')=='plot_survival','Wrong survival mapping')
                    require('constplot' in options,'Empirical survival must retain its step display')
                    kind='survival';endpoint_values=lambda row:[row['plot_survival']]
                else:
                    power=role[-1]
                    require(match.group(1)=='moments' and columns.get('y')=='lo'+power,'Wrong moment mapping')
                    if columns.get('y error')=='err'+power:
                        require('ydir=plus' in options,'Computational range must start at the lower endpoint')
                        kind='computational'
                        endpoint_values=lambda row,p=power:[row['lo'+p],row['hi'+p]]
                    else:
                        require(columns.get('y error plus')=='outerplus'+power
                                and columns.get('y error minus')=='outerminus'+power and 'ydir=both' in options,
                                'Wrong approximate sampling-envelope columns')
                        kind='sampling'
                        endpoint_values=lambda row,p=power:[row['outerlo'+p],row['outerhi'+p]]
                if kind!='sampling': require('mark='+MARKS[mixture]+',' in options,'Changed mixture marker identity')
                key=(mixture,kind);observed[key]=observed.get(key,0)+1
                for row in rows:
                    for y in endpoint_values(row):
                        if math.isnan(y): continue
                        require(bounds['xmin']<=row['x']<=bounds['xmax']
                                and bounds['ymin']-1e-12<=y<=bounds['ymax']+1e-12,
                                'A retained plotted endpoint is clipped')
                        endpoint_count+=1
            kinds=['failure'] if role in ('g','g_prime') else ['survival'] if role=='survival' else ['computational','sampling']
            require(observed=={(m,k):1 for m in TYPES for k in kinds},'Missing or duplicate mixture series')
            panels.append({'role':role,'axis_bounds':bounds,'table_commands':len(commands),
                           'visible_endpoint_values_checked':endpoint_count,'all_required_series_present':True})
    expected_roles= ['g','g_prime'] if path.name=='fig2_certificate_failure.tex' else ['survival','moment1','moment2'] if path.name=='fig3_residual_and_moments.tex' else ['g','g_prime','survival','moment1','moment2']
    require([p['role'] for p in panels]==expected_roles,'Wrong scientific-panel inventory/order')
    expected_groups=[[2,1]] if len(panels)==2 else [[3,1]] if len(panels)==3 else [[2,1],[3,1]]
    require(group_sizes==expected_groups,'The optional overview must retain two panels above three')
    return {'source_sha256':digest(path),'group_sizes':group_sizes,'panels':panels}


def check_zero_policy(figures):
    source=(figures/'build_bernoulli_figures.py').read_text()
    compact=re.sub(r'\s+','',source)
    tokens=['p=pifcountelse"nan"','upper=hi','mark=triangle','rotate=180','fill=white','x=x,y=upper','\\triangledown']
    for token in tokens:
        require(token in compact,'Reusable zero-failure policy missing: '+token)
    return {'status':'PASS','method':'Source-policy check; all retained failure counts are positive',
            'behavior':'Zero count retains true rate zero, uses nan on the log-frequency curve, and marks the exact CP upper endpoint with an open downward triangle.'}


def check_presentation_metadata(figures,result):
    path=figures/'condition_moment_presentation.json'
    metadata=json.loads(path.read_text())
    require(metadata['network_sizes']==SIZES and metadata['network_size_display_offsets']==OFFSETS,
            'Presentation size grid or offsets disagree')
    require(metadata['mixture_colors']==COLORS and metadata['mixture_markers']==MARKS,
            'Presentation mixture identities disagree')
    require(metadata['active_figures']==['fig2_certificate_failure','fig3_residual_and_moments'],
            'Unexpected active figure inventory')
    require(metadata['optional_consolidation']=='fig_condition_moment_overview',
            'Unexpected optional consolidation identity')
    require((figures/'fig_condition_moment_overview.tex').is_file(),
            'Declared optional consolidation source is missing')
    require(metadata['failure_axes']=={'scale':'log-log','x':[70,36000],'y':[.0000015,.7]},
            'Failure-axis metadata disagrees')
    require(metadata['survival_axes']=={'scale':'linear-log','x':[0,200],'y':[.00007,1.2]},
            'Survival-axis metadata disagrees')
    for power,expected in [('1',[.8275195516154046,1.2107929637406762]),
                           ('2',[.5315788186407935,3.289914264413321])]:
        for actual,reference in zip(metadata['moment_axes'][power],expected):
            close(actual,reference,'Actual moment-axis metadata')
        for actual,reference in zip(metadata['moment_outer_envelope_extrema'][power],
                                    result['normalized_moment_outer_endpoint_extrema'][int(power)]):
            close(actual,reference,'Moment-envelope metadata')
    require(metadata['retained_plot_hashes']=={name:value for name,value in CSV_SHA256.items()
            if name.startswith('bernoulli_')},'Presentation/data integrity hashes disagree')
    for name,expected in metadata['tex_hashes'].items():
        require(digest(figures/(name+'.tex'))==expected,'Presentation/source hash mismatch: '+name)
    require(set(metadata['tex_hashes'])=={'fig2_certificate_failure','fig3_residual_and_moments',
                                        'fig_condition_moment_overview'},'Source hash inventory mismatch')
    return {'status':'PASS','sha256':digest(path),'axes_and_envelope_extrema_distinguished':True,
            'retained_plot_hashes_verified':len(metadata['retained_plot_hashes']),
            'tex_hashes_verified':len(metadata['tex_hashes'])}


def check_all(figures):
    result=check_data(figures.parent)
    result['figure_source_checks']={}
    for name in ['fig2_certificate_failure.tex','fig3_residual_and_moments.tex','fig_condition_moment_overview.tex']:
        path=figures/name
        if path.exists(): result['figure_source_checks'][name]=inspect_figure(path)
        else:
            require(name=='fig_condition_moment_overview.tex','A required separate figure is missing')
            result['overview_check']='NOT INCLUDED: separate figures checked'
    result['zero_failure_policy']=check_zero_policy(figures)
    result['presentation_metadata_check']=check_presentation_metadata(figures,result)
    result['source_scope']='Included aggregate tables and actual plotted TeX; no network experiment or approximate-SE replay'
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='JSON destination; defaults to source-root qa/')
    args=parser.parse_args()
    figures=Path(__file__).resolve().parent
    output=args.output or figures.parent/'qa/condition_moment_verification.json'
    result=check_all(figures)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['status','failure_rows_checked','cp_endpoints_recomputed',
         'cp_recalculation','moment_rows_checked','normalized_moment_ranges_checked','survival_rows_checked',
         'zero_survival_rows_retained','csv_files_checked']},indent=2))
