"""Shared isolated computational-figure checks; no simulations are run."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
from reproduce import RunnerError, _output_destination, assertion_policy, verify_repository

ROOT = Path(__file__).resolve().parents[1]
ACTIVE = (
    ('fig1_cost_distribution', 'Cost probabilities and Poisson comparison'),
    ('fig2_certificate_failure', 'Copy-condition failure frequencies'),
    ('fig3_residual_and_moments', 'Residual survival and normalized moments'),
    ('fig5_larger_network_support_grouped', 'Support bounds and grouped original-state orbit'),
    ('fig6_cd4_target_mechanisms', 'Complete target costs and isolated component maps'),
    ('fig4_biological_certificates', 'Two exact synchronous mechanisms'),
)
COMPANIONS = (
    'fig1_cost_differences',
    'fig5_larger_network_support', 'fig_condition_moment_overview',
)
GENERATED = (
    'fig1_cost_distribution.tex', 'fig1_cost_differences.tex',
    'fig2_certificate_failure.tex', 'fig3_residual_and_moments.tex',
    'fig_condition_moment_overview.tex', 'cost_distribution_presentation.json',
    'condition_moment_presentation.json',
)
VALIDATORS = (
    'verify_cost_distribution.py', 'verify_condition_moment_data.py',
    'verify_regulatory_figures.py', 'verify_grouped_support.py',
)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def snapshot(path):
    return {p.relative_to(path).as_posix(): sha(p) for p in sorted(path.rglob('*')) if p.is_file()}

def require(condition, explanation):
    if not condition:
        raise RunnerError(explanation)

def json_write(path, record):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')

def check_repository(expected=None):
    require(assertion_policy()['enabled'], 'Run without optimized Python; assertions must remain enabled.')
    report = verify_repository(ROOT, expected_manifest_sha256=expected)
    require(report['ok'], 'Repository integrity check failed: '+'; '.join(report['errors']))
    return report

def run(argv, cwd, log, timeout=180):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', SOURCE_DATE_EPOCH='1792022400', FORCE_SOURCE_DATE='1')
    result = subprocess.run([str(x) for x in argv], cwd=cwd, env=env, capture_output=True, text=True, errors='replace', timeout=timeout)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(result.stdout+result.stderr)
    require(result.returncode == 0, 'Command failed; see '+str(log))
    return result.stdout

def check_inventory(figures):
    record=json.loads((figures/'FIGURE_INVENTORY.json').read_text())
    require(record['schema']==1, 'Unsupported figure-inventory schema')
    require(isinstance(record.get('files'),dict) and bool(record['files']),
            'Figure inventory is empty or invalid')
    for relative, entry in record['files'].items():
        require(not Path(relative).is_absolute() and '..' not in Path(relative).parts,
                'Unsafe figure inventory path: '+relative)
        path=figures/relative
        require(path.is_file() and path.stat().st_size==entry['bytes'] and sha(path)==entry['sha256'],
                'Figure inventory mismatch: '+relative)
    metadata=json.loads((figures/'bernoulli_plot_metadata.json').read_text())
    source_directory=Path(metadata['source_directory'])
    require(not source_directory.is_absolute() and '..' not in source_directory.parts,
            'Unsafe statistical source directory')
    require(isinstance(metadata.get('input_hashes'),dict) and bool(metadata['input_hashes']),
            'Statistical input inventory is empty or invalid')
    for name,expected in metadata['input_hashes'].items():
        require(not Path(name).is_absolute() and '..' not in Path(name).parts,
                'Unsafe statistical input path: '+name)
        require(sha(ROOT/source_directory/name)==expected,
                'Statistical export input differs: '+name)
    return {'files':len(record['files']),'statistical_inputs_checked':len(metadata['input_hashes']),
            'scope':'Current figure assets and numerical input digests.'}

def verify_copy(figures, reports):
    """Check copied retained inputs and regenerate only copied numerical sources."""
    reports.mkdir(parents=True, exist_ok=True)
    original=snapshot(figures)
    inventory=check_inventory(figures)
    checks={}
    for relative in VALIDATORS:
        key=Path(relative).stem
        output=reports/(key+'.json')
        run([sys.executable, figures/relative, '--output', output], figures,
            reports/(key+'.log'))
        result=json.loads(output.read_text())
        # Individual validators may return structured scientific checks without
        # a global status; their nonzero exit is still mandatory on failure.
        require(result.get('status','PASS')=='PASS', 'Verification reported failure: '+key)
        checks[key]=result
    require(snapshot(figures)==original, 'A verifier altered copied input assets')
    run([sys.executable, figures/'build_bernoulli_figures.py', '--sources-only'], figures,
        reports/'generator_agreement.log')
    regenerated=snapshot(figures)
    mismatches=[name for name in GENERATED if regenerated.get(name)!=original.get(name)]
    require(not mismatches, 'Generator-to-source agreement failed: '+', '.join(mismatches))
    require(regenerated==original, 'Source regeneration changed an unexpected copied input')
    return {'status':'PASS','inventory':inventory,'checks':checks,
            'generator_agreement':{'tex_sources':5,'presentation_records':2,'all_input_bytes_unchanged':True},
            'scope':'Retained exports, native plot mappings and small schematic/component maps; no scientific workflow or raw-network replay.'}

def verify_current():
    before=check_repository()
    try:
        with tempfile.TemporaryDirectory(prefix='boolean_pinning_figures_') as name:
            temporary=Path(name)
            figures=temporary/'figures'
            shutil.copytree(ROOT/'figures',figures)
            report=verify_copy(figures,temporary/'checks')
    finally:
        after=check_repository(before['manifest_sha256'])
    report['repository_integrity_before']=before
    report['repository_integrity_after']=after
    return report

def compare_pdf(reference, rebuilt, scratch):
    scratch.mkdir(parents=True,exist_ok=True)
    info=[];text=[];pixels=[]
    for label,path in [('reference',reference),('rebuilt',rebuilt)]:
        listing=run(['pdfinfo',path],scratch,scratch/(label+'_pdfinfo.log'))
        pages=re.search(r'^Pages:\s+(\d+)',listing,re.M)
        require(pages is not None and int(pages.group(1))==1, 'Expected a complete one-page native figure: '+str(path))
        require(not re.search(r'^Author:[ \t]*\S',listing,re.M),'Unexpected personal author metadata in '+str(path))
        textpath=scratch/(label+'.txt')
        run(['pdftotext','-layout',path,textpath],scratch,scratch/(label+'_text.log'))
        prefix=scratch/label
        run(['pdftoppm','-singlefile','-r','72',path,prefix],scratch,scratch/(label+'_render.log'))
        text.append(textpath.read_bytes());pixels.append(prefix.with_suffix('.ppm').read_bytes())
        info.append(sha(path))
    require(text[0]==text[1], 'Rebuilt figure text differs: '+reference.name)
    require(pixels[0]==pixels[1], 'Rebuilt figure pixels/dimensions differ: '+reference.name)
    return {'reference_sha256':info[0],'rebuilt_sha256':info[1],'byte_identical':info[0]==info[1],
            'text_identical':True,'dimensions_and_72dpi_pixels_identical':True}

def build_current(output):
    destination=_output_destination(ROOT,Path(output))
    before=check_repository()
    destination.mkdir(parents=True,exist_ok=False)
    report={'status':'FAIL','scope':'Current computational figures only; no simulations or model census.',
            'repository_integrity_before':before,'source_date_epoch':1792022400,'latex_passes':3}
    try:
        figures=destination/'workspace/figures'
        shutil.copytree(ROOT/'figures',figures)
        report['verification']=verify_copy(figures,destination/'checks')
        inputs=snapshot(figures)
        for command in ['pdflatex','pdfinfo','pdftotext','pdftoppm']:
            require(shutil.which(command) is not None, 'Required figure tool is unavailable: '+command)
        outputs={}
        for relative in [name for name,_ in ACTIVE]+list(COMPANIONS):
            source=figures/(relative+'.tex')
            require(source.is_file(),'Missing native figure source: '+relative)
            for pass_number in range(1,4):
                run(['pdflatex','-interaction=nonstopmode','-halt-on-error',source.name],source.parent,
                    destination/'logs'/(source.stem+f'_pass{pass_number}.log'))
            log=source.with_suffix('.log').read_text(errors='replace')
            require(not any(token in log for token in ['Overfull','undefined references','undefined citations','LaTeX Error:']),
                    'Unresolved or overflowing native figure: '+relative)
            built=source.with_suffix('.pdf')
            comparison=compare_pdf(ROOT/'figures'/(relative+'.pdf'),built,
                                   destination/'comparison'/source.stem)
            target=destination/'figures'/(relative+'.pdf')
            target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(built,target)
            outputs[relative]=comparison
        # The copied PDFs and compiler auxiliaries may change; every shipped
        # non-PDF source/data/metadata input must retain its original bytes.
        current=snapshot(figures)
        changed=[p for p,h in inputs.items() if not p.endswith('.pdf') and current.get(p)!=h]
        require(not changed,'Build altered copied source or data: '+', '.join(changed))
        report.update(status='PASS',outputs=outputs,active_figures=len(ACTIVE),companion_figures=len(COMPANIONS),
                      copied_non_pdf_inputs_preserved=True)
    except Exception as exc:
        report['error']=str(exc)
        raise
    finally:
        try:
            report['repository_integrity_after']=check_repository(before['manifest_sha256'])
        except Exception as exc:
            report['status']='FAIL';report['integrity_error']=str(exc)
            json_write(destination/'build_summary.json',report)
            raise
        json_write(destination/'build_summary.json',report)
    return report
