"""d02: model sizes -- atoms, electrons and basis-function counts recomputed with PySCF from the case inputs.

Run with the DFT env: /root/cqc-sfa-tools/env/bin/python d02_models.py
Spherical (5d) functions as in the production runs; 6-311++G(d,p) unless the case names another basis.
Also counts the functions actually used after PySCF's linear-dependency removal, as logged by the SCF runs
(n_basis in fmo.json / opt_summary.json) -- read, not recomputed.
"""
import json
from pyscf import gto
from common import *

CASES = ['sfa_pair_pcm_d3', 'sfa_pair_gas_d3', 'sfa_pair_pcm_b3lyp', 'sfa_pair_pcm_d3_tz', 'cation_pcm_d3',
         'cation_gas_d3', 'h2po4_pcm_d3', 'h2po4_gas_d3', 'h3po4_pcm_d3', 'h3po4_gas_d3', 'sulfanilamide_pcm_d3',
         'sulfanilamide_gas_d3', 'sulfanilamide_pcm_d3_tz', 'sulfonamidate_pcm_d3', 'urea_gas_d3']
out = {}
for c in CASES:
    case = json.loads((REPO / f'molecules/sfa/cases/{c}.json').read_text())
    basis = case.get('basis', '6-311++g(d,p)')
    mol = gto.M(atom=str(REPO / case['xyz']), basis=basis, charge=case['charge'], spin=case['spin'],
                cart=False, verbose=0)
    syms = [mol.atom_pure_symbol(i) for i in range(mol.natm)]
    formula = {}
    for s in syms:
        formula[s] = formula.get(s, 0) + 1
    out[c] = dict(xyz=case['xyz'], basis=basis, charge=case['charge'], multiplicity=case['spin'] + 1,
                  n_atoms=mol.natm, n_heavy=sum(1 for s in syms if s != 'H'), n_electrons=mol.nelectron,
                  nao_spherical=mol.nao_nr(), formula=formula, solvent=case.get('solvent'), disp=case.get('disp'))
# the (pair)2 cluster planned as M8 was never built; its size is simply twice the pair
out['_note'] = 'M8 (pair)2 cluster of PLAN.md was not built or run; no result exists for it.'
dump('d02_models.json', out)
