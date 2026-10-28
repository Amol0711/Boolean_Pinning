"""Check current protocol, seed, validation and production-source identities."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

CURRENT_FILES = {
    'protocol_sha256': 'Protocol_J1_Replication.json',
    'seed_manifest_sha256': 'seed_manifest.csv',
    'independent_preflight_sha256': 'qa/independent_preflight.json',
    'production_preflight_sha256': 'qa/production_preflight.json',
    'independent_oracle_sha256': 'verification/independent_network_oracle.py',
}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_representation(root):
    """Verify current file identities against both integrity inventories."""
    root = Path(root)
    current = json.loads((root/'integrity_record.json').read_text())
    if current.get('schema') != 1:
        raise ValueError('Unsupported integrity record')
    files = {name: current[field] for field, name in CURRENT_FILES.items()}
    files.update(current['source_sha256'])
    files['integrity_record.json'] = sha(root/'integrity_record.json')
    files['representation.py'] = sha(root/'representation.py')
    inventory = json.loads((root/'SHA256_MANIFEST.json').read_text())['files']
    for name, expected in files.items():
        if sha(root/name) != expected or inventory.get(name) != expected:
            raise ValueError('File integrity mismatch: '+name)
    return current, current

def protocol_reference(record):
    """Describe the protocol checked for a computation."""
    return {'protocol_sha256': record['protocol_sha256'],
            'seed_manifest_sha256': record['seed_manifest_sha256'],
            'scope': 'Current protocol and fixed random-stream identities.'}
