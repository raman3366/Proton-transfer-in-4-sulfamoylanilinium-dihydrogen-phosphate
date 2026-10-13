#!/usr/bin/env python3
"""
ADMET / drug-likeness sanity check for the SI (user decision 27 Sep 2026: SI only, as a sanity check against
sulfanilamide's known clinical and experimental profile).  Everything is computed locally; no web form is used.

Usage (repository root):
    nice -n 10 /root/cqc-sfa-tools/admet/bin/python paper/si/compute/admet.py [--tdc-dir DIR]
Writes molecules/sfa/results/sulfanilamide_pcm_d3/admet/{admet.json, admet.csv} and the same pair for the cation
under molecules/sfa/results/cation_pcm_d3/admet/ (reference only: the cation does not exist at pH 7.4, and the
ADMET-AI models were trained on neutral drug SMILES, so its predictions are out of domain and flagged as such).

What is computed
  1. Identity: the SMILES used here is checked against the DFT structure of the same case: bonds perceived from
     results/<case>/opt.xyz with RDKit rdDetermineBonds at the case charge must give the same InChIKey.
  2. RDKit descriptors: MolWt, exact mass, heavy/total atoms, Wildman-Crippen logP and MR, TPSA (Ertl; N and O
     only, as in the original method, and with S and P), H-bond donors/acceptors (Lipinski's NH+OH and N+O counts,
     and RDKit's pharmacophore-based counts), rotatable bonds (RDKit default, strict), rings, Fsp3, formal charge,
     QED.  Rule filters with their published thresholds: Lipinski (MW <= 500, logP <= 5, HBD <= 5, HBA <= 10),
     Veber (rotatable bonds <= 10, TPSA <= 140), Ghose (logP -0.4..5.6, MR 40..130, MW 160..480, 20..70 atoms
     including H), Egan (logP <= 5.88, TPSA <= 131.6).  Egan's and Ghose's own logP were AlogP98 / ALOGP;
     Wildman-Crippen logP is used for all of them and this is stated.
  3. ADMET-AI 1.4.0 (Swanson et al. 2024, Chemprop-RDKit ensembles on the TDC ADMET group), with the percentile
     of each value among the 2579 approved DrugBank drugs that ADMET-AI ships as its reference set.
  4. Training-set membership (optional, --tdc-dir with the TDC .tab files downloaded from Harvard Dataverse):
     for every TDC dataset behind an ADMET-AI endpoint, whether sulfanilamide is in it and with what label.  A
     prediction for a training molecule is a recall of its label, not an independent prediction.
"""
import argparse, datetime, hashlib, json, os, platform, sys
from pathlib import Path

os.environ.setdefault('OMP_NUM_THREADS', '2')
REPO = Path(__file__).resolve().parents[3]
RES = REPO / 'molecules/sfa/results'

ap = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
ap.add_argument('--tdc-dir', help='directory with TDC <dataset>.tab files (Harvard Dataverse downloads)')
ap.add_argument('--tdc-ids', help='JSON {dataset: [dataverse file id, type]} used to download them')
args = ap.parse_args()

import numpy as np                                                    # noqa: E402
import torch                                                          # noqa: E402
torch.set_num_threads(2)
import rdkit                                                          # noqa: E402
from rdkit import Chem, RDLogger                                      # noqa: E402
from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors, rdDetermineBonds  # noqa: E402
RDLogger.DisableLog('rdApp.*')
import admet_ai, chemprop                                             # noqa: E402
from admet_ai import ADMETModel                                       # noqa: E402

SPECIES = {
    'sulfanilamide_pcm_d3': {'name': 'sulfanilamide (neutral)', 'smiles': 'Nc1ccc(cc1)S(N)(=O)=O', 'charge': 0,
                             'role': 'primary: the species present at pH 7.4 (anilinium pKa about 2, sulfonamide pKa about 10.4)'},
    'cation_pcm_d3': {'name': '4-sulfamoylanilinium cation', 'smiles': '[NH3+]c1ccc(cc1)S(N)(=O)=O', 'charge': 1,
                      'role': 'reference only: not present at pH 7.4; ADMET-AI endpoints are out of their training domain'},
}


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def identity(case, sp):
    mol = Chem.MolFromSmiles(sp['smiles'])
    xyz = RES / case / 'opt.xyz'
    m3 = Chem.MolFromXYZFile(str(xyz))
    rdDetermineBonds.DetermineBonds(m3, charge=sp['charge'])
    m3 = Chem.RemoveHs(m3)
    ik_smi, ik_xyz = Chem.MolToInchiKey(mol), Chem.MolToInchiKey(m3)
    return mol, {'smiles_input': sp['smiles'], 'smiles_canonical': Chem.MolToSmiles(mol),
                 'inchikey': ik_smi, 'inchi': Chem.MolToInchi(mol),
                 'dft_geometry': str(xyz.relative_to(REPO)), 'dft_geometry_sha256': sha256(xyz),
                 'smiles_from_dft_geometry': Chem.MolToSmiles(m3), 'inchikey_from_dft_geometry': ik_xyz,
                 'identity_ok': ik_smi == ik_xyz}


def descriptors(mol):
    mh = Chem.AddHs(mol)
    d = {
        'MW': Descriptors.MolWt(mol), 'exact_mass': Descriptors.ExactMolWt(mol),
        'formula': rdMolDescriptors.CalcMolFormula(mol), 'formal_charge': Chem.GetFormalCharge(mol),
        'heavy_atoms': mol.GetNumHeavyAtoms(), 'atoms_incl_H': mh.GetNumAtoms(),
        'cLogP_Crippen': Crippen.MolLogP(mol), 'MR_Crippen': Crippen.MolMR(mol),
        'TPSA_NO': rdMolDescriptors.CalcTPSA(mol), 'TPSA_NOSP': rdMolDescriptors.CalcTPSA(mol, includeSandP=True),
        'HBD_Lipinski_NHOH': rdMolDescriptors.CalcNumLipinskiHBD(mol), 'HBA_Lipinski_NO': rdMolDescriptors.CalcNumLipinskiHBA(mol),
        'HBD_RDKit': rdMolDescriptors.CalcNumHBD(mol), 'HBA_RDKit': rdMolDescriptors.CalcNumHBA(mol),
        'rotatable_bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
        'rings': rdMolDescriptors.CalcNumRings(mol), 'aromatic_rings': rdMolDescriptors.CalcNumAromaticRings(mol),
        'Fsp3': rdMolDescriptors.CalcFractionCSP3(mol), 'QED': QED.qed(mol),
    }
    rules = {
        'Lipinski': {'MW <= 500': d['MW'] <= 500, 'cLogP <= 5': d['cLogP_Crippen'] <= 5,
                     'HBD (NH+OH) <= 5': d['HBD_Lipinski_NHOH'] <= 5, 'HBA (N+O) <= 10': d['HBA_Lipinski_NO'] <= 10},
        'Veber': {'rotatable bonds <= 10': d['rotatable_bonds'] <= 10, 'TPSA <= 140': d['TPSA_NO'] <= 140},
        'Ghose': {'-0.4 <= cLogP <= 5.6': -0.4 <= d['cLogP_Crippen'] <= 5.6, '40 <= MR <= 130': 40 <= d['MR_Crippen'] <= 130,
                  '160 <= MW <= 480': 160 <= d['MW'] <= 480, '20 <= atoms (incl. H) <= 70': 20 <= d['atoms_incl_H'] <= 70},
        'Egan': {'cLogP <= 5.88': d['cLogP_Crippen'] <= 5.88, 'TPSA <= 131.6': d['TPSA_NO'] <= 131.6},
    }
    summary = {k: {'violations': [c for c, ok in v.items() if not ok], 'n_violations': sum(not ok for ok in v.values()),
                   'passes': all(v.values()) if k != 'Lipinski' else sum(not ok for ok in v.values()) <= 1}
               for k, v in rules.items()}
    return d, rules, summary


model = ADMETModel(num_workers=0)
admet_info = {r['id']: r for r in __import__('csv').DictReader(open(Path(admet_ai.__file__).parent / 'resources/data/admet.csv'))}
drugbank = model.drugbank

tdc = None
if args.tdc_dir:
    import csv, glob
    ids = json.loads(Path(args.tdc_ids).read_text()) if args.tdc_ids else {}
    tdc = {}
    target = Chem.MolToInchiKey(Chem.MolFromSmiles(SPECIES['sulfanilamide_pcm_d3']['smiles']))
    for f in sorted(glob.glob(os.path.join(args.tdc_dir, '*.tab'))):
        name = Path(f).stem
        rows = list(csv.reader(open(f), delimiter='\t'))
        hdr = rows[0]
        smi_col = [i for i, h in enumerate(hdr) if h in ('Drug', 'X')]
        hits = []
        for r in rows[1:]:
            for i in smi_col:
                m = Chem.MolFromSmiles(r[i])
                if m is not None and Chem.MolToInchiKey(m) == target:
                    hits.append(dict(zip(hdr, r)))
        tdc[name] = {'n': len(rows) - 1, 'sha256': sha256(f), 'dataverse_file_id': (ids.get(name) or [None])[0],
                     'sulfanilamide_in_set': bool(hits), 'rows': hits}

stamp = datetime.datetime.now().isoformat(timespec='seconds')
versions = {'python': platform.python_version(), 'rdkit': rdkit.__version__, 'admet_ai': admet_ai.__version__
            if hasattr(admet_ai, '__version__') else '1.4.0', 'chemprop': chemprop.__version__, 'torch': torch.__version__,
            'numpy': np.__version__, 'machine': platform.machine(), 'host': platform.node()}
TDC_OF = {k.lower(): k for k in admet_info}          # tdc file stem (lower case) -> ADMET-AI endpoint id

for case, sp in SPECIES.items():
    mol, ident = identity(case, sp)
    if not ident['identity_ok']:
        raise SystemExit(f'[admet] {case}: SMILES and DFT geometry disagree: {ident}')
    d, rules, summ = descriptors(mol)
    pred = model.predict(smiles=ident['smiles_canonical'])
    pred = {k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in pred.items()}
    endpoints = {}
    for k, info in admet_info.items():
        e = {'name': info['name'], 'category': info['category'], 'task': info['task_type'], 'units': info['units'],
             'value': pred.get(k), 'drugbank_percentile': pred.get(f'{k}_drugbank_approved_percentile')}
        # membership refers to neutral sulfanilamide; the cation (a different InChIKey) is in none of the sets
        if tdc is not None and info['category'] != 'Physicochemical' and case == 'sulfanilamide_pcm_d3':
            stem = k.lower() if k.lower() in tdc else ('tox21' if k.startswith(('NR-', 'SR-')) else None)
            if stem:
                t = tdc[stem]
                lab = None
                if t['rows']:
                    row = t['rows'][0]
                    lab = row.get(k) if stem == 'tox21' else row.get('Y')
                e['training_set'] = {'dataset': stem, 'sulfanilamide_in_set': t['sulfanilamide_in_set'],
                                     'label': (float(lab) if lab not in (None, '') else None)}
        endpoints[k] = e
    stored = None
    if drugbank is not None and case == 'sulfanilamide_pcm_d3':
        row = drugbank[drugbank['smiles'].map(lambda s: Chem.MolToInchiKey(Chem.MolFromSmiles(s)) == ident['inchikey'])]
        if len(row):
            r0 = row.iloc[0]
            diffs = {k: abs(float(r0[k]) - endpoints[k]['value']) for k in admet_info if k in r0 and endpoints[k]['value'] is not None}
            stored = {'drugbank_id': r0.get('id'), 'max_abs_diff_vs_shipped_prediction': max(diffs.values()),
                      'n_compared': len(diffs)}
    out = {'provenance': {'script': 'paper/si/compute/admet.py', 'argv': sys.argv[1:], 'finished': stamp, 'versions': versions,
                          'admet_ai_reference_set': 'DrugBank approved drugs shipped with ADMET-AI (drugbank_approved.csv), '
                                                    f'{len(drugbank) if drugbank is not None else 0} molecules',
                          'threads': 2},
           'species': sp['name'], 'role': sp['role'], 'identity': ident,
           'rdkit_descriptors': d, 'rules': {k: {c: bool(ok) for c, ok in v.items()} for k, v in rules.items()},
           'rule_summary': summ, 'admet_ai': endpoints, 'admet_ai_vs_shipped_drugbank_row': stored,
           'notes': ['cLogP is Wildman-Crippen (RDKit); Ghose and Egan defined their filters with ALOGP/AlogP98.',
                     'TPSA_NO is Ertl TPSA over N and O (the original definition, used for Veber/Egan); TPSA_NOSP adds S and P.',
                     'ADMET-AI classification outputs are model probabilities (0-1), not measured outcomes.',
                     'A prediction for a molecule that is in the training set reproduces its label; see training_set.']}
    if tdc is not None and case == 'sulfanilamide_pcm_d3':
        out['tdc_training_sets'] = {k: {kk: vv for kk, vv in v.items() if kk != 'rows'} | {'rows': v['rows']} for k, v in tdc.items()}
    od = RES / case / 'admet'
    od.mkdir(parents=True, exist_ok=True)
    (od / 'admet.json').write_text(json.dumps(out, indent=2, default=str))
    lines = ['group,key,value,unit,drugbank_percentile,in_training_set,training_label']
    for k, v in d.items():
        lines.append(f'rdkit,{k},{v},,,,')
    for k, e in endpoints.items():
        ts = e.get('training_set') or {}
        lines.append(f"admet_ai,{k},{e['value']},{e['units']},{e['drugbank_percentile']},{ts.get('sulfanilamide_in_set', '')},{ts.get('label', '')}")
    (od / 'admet.csv').write_text('\n'.join(lines) + '\n')
    print(f"[admet] {case}: identity {ident['identity_ok']} ({ident['smiles_canonical']}); MW {d['MW']:.2f}, cLogP "
          f"{d['cLogP_Crippen']:.3f}, TPSA {d['TPSA_NO']:.2f}; rules " +
          ', '.join(f"{k} {v['n_violations']} viol." for k, v in summ.items()) +
          (f"; vs shipped DrugBank row max |diff| {stored['max_abs_diff_vs_shipped_prediction']:.2e}" if stored else ''))
