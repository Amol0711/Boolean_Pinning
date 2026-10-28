#!/usr/bin/env python3
"""Independent source-corner audit; does not import extraction/engine modules."""
from pathlib import Path
import ast,hashlib,itertools,json
from fractions import Fraction
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def index(node,array):
    if isinstance(node,ast.Subscript) and isinstance(node.value,ast.Name) and node.value.id==array and isinstance(node.slice,ast.Constant) and type(node.slice.value) is int:return node.slice.value
    return None

def calculate(node,x):
    if isinstance(node,ast.Constant) and type(node.value) in (int,float):return Fraction(str(node.value))
    i=index(node,'S')
    if i is not None:return Fraction(x[i])
    if isinstance(node,ast.BinOp):
        a,b=calculate(node.left,x),calculate(node.right,x)
        if isinstance(node.op,ast.Add):return a+b
        if isinstance(node.op,ast.Sub):return a-b
        if isinstance(node.op,ast.Mult):return a*b
    if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return -calculate(node.operand,x)
    raise ValueError(ast.dump(node))

def bnet_value(node,x):
    if isinstance(node,ast.Name):return x[int(node.id[1:])-1]
    if isinstance(node,ast.Constant) and node.value in (0,1):return node.value
    if isinstance(node,ast.BoolOp):
        v=[bnet_value(t,x) for t in node.values]
        if isinstance(node.op,ast.And):return int(all(v))
        if isinstance(node.op,ast.Or):return int(any(v))
    if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.Not):return 1-bnet_value(node.operand,x)
    raise ValueError(ast.dump(node))

records=[]
for n,stem in [(17,'drosophila_coded_17'),(29,'tlgl_coded_29')]:
    source=HERE/'source_files'/f'{n}nodeQL.py';docpath=ROOT/'models'/f'{stem}.json';bnetpath=docpath.with_suffix('.bnet')
    tree=ast.parse(source.read_bytes().decode());doc=json.loads(docpath.read_text())
    fn=next(t for t in tree.body if isinstance(t,ast.FunctionDef) and t.name=='calculate_H')
    expressions={index(t.targets[0],'new_S'):t.value for t in ast.walk(fn) if isinstance(t,ast.Assign) and len(t.targets)==1 and index(t.targets[0],'new_S') is not None}
    assert set(expressions)==set(range(n))
    targetnode=next(t.value for t in tree.body if isinstance(t,ast.Assign) and any(isinstance(x,ast.Name) and x.id=='attractor' for x in t.targets))
    target=[int(t) for t in ast.literal_eval(targetnode)]
    assert doc['target']==target
    assert all(calculate(expressions[i],target)==target[i] for i in range(n))
    bro=[s for s in bnetpath.read_text().splitlines() if s.strip() and not s.startswith('#') and not s.startswith('targets,')]
    assert len(bro)==n
    source_count=0;union_count=0
    for i,r in enumerate(doc['rules']):
        assert r['node_zero_based']==i and r['node_one_based']==i+1
        expr=expressions[i];scope=sorted({index(q,'S') for q in ast.walk(expr) if index(q,'S') is not None})
        parents=r['parents_zero_based'];assert len(parents)==len(set(parents));assert len(r['physical_truth_table'])==2**len(parents)
        label,formula=bro[i].split(',',1);assert label.strip()==f'n{i+1:02d}'
        bt=ast.parse(formula.replace('!',' not ').replace('&',' and ').replace('|',' or ').strip(),mode='eval').body
        bparents={int(q.id[1:])-1 for q in ast.walk(bt) if isinstance(q,ast.Name)}
        union=sorted(set(scope)|set(parents)|bparents)
        for vals in itertools.product((0,1),repeat=len(union)):
            x=[0]*n
            for j,v in zip(union,vals):x[j]=v
            original=calculate(expr,x)
            exported=r['physical_truth_table'][sum(x[j]<<k for k,j in enumerate(parents))]
            bnet=bnet_value(bt,x)
            assert original in (0,1) and original==exported==bnet,(n,i,x)
            union_count+=1
        source_count+=2**len(scope)
    records.append({'model_id':stem,'n':n,'source_sha256':sha(source),'json_sha256':sha(docpath),'bnet_sha256':sha(bnetpath),'source_local_corner_count':source_count,'union_corners_checked':union_count,'all_original_source_JSON_BNET_corners_identical':True,'target_equals_original_source_literal':True,'source_target_fixed':True,'control_node_primary_labels':{'1':'SLP','2':'nWG','3':'nhhnHH (source alias nHH)'} if n==17 else {'1':'IL15','9':'PDGF','15':'SPHK','16':'S1P'}})
result={'status':'PASS','method':'Original preserved frozen source expressions parsed directly; exact Fraction AST interpreter independent of package imports; physical JSON table lookup and independent BNET parser compared on every union-scope corner; no prior registry used','models':records}
(HERE/'independent_export_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
