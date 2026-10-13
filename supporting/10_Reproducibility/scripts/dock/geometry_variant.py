"""
Write a diagnostic ligand geometry that differs from a reference geometry in exactly one internal
coordinate, so that the effect of that coordinate on a docking result can be measured.

    python dock/geometry_variant.py --species sfa_neutral \
        --xyz molecules/sfa/results/sulfanilamide_pcm_d3/opt.xyz \
        --set-angle C4 S1 N2 102.0 --out dock/results/diagnostics/sfa_neutral_csn102.xyz

Why this exists (production 2026-09-28).  Docking the DFT geometry instead of the ETKDG/MMFF94s
conformer moved the AutoDock4Zn pose of sulfanilamide by ~3 Å in 3HS4 and removed the bimodality in
6RL9, with the zinc anchor unchanged in both.  The two inputs differ most in the C4-S1-N2 angle
(MMFF94s 102.0 deg, B3LYP-D3/PCM 109.5 deg, crystal structures 106.9-108.8 deg), and Vina keeps bond
angles rigid, so that angle is the obvious candidate.  Changing it on its own, in both directions,
turns the observation into a test.

Atom names are the crystal labels of this compound (C1..C6, N1, N2, O1, O2, S1), which the DFT xyz
files follow; `--index` takes 1-based positions instead for any other file.  The written xyz keeps the
input atom order, so `prep_ligands.from_dft_xyz` reads it exactly as it reads the DFT file.
"""
import argparse
import pathlib
import sys

from rdkit import Chem
from rdkit.Chem import rdMolTransforms as T

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from prep_ligands import SPECIES, from_dft_xyz  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
# the label -> 0-based index convention of the sulfanilamide DFT xyz files (checked in
# production_analysis.dft_mol, which asserts the bonding these labels imply)
LABELS = dict(C1=0, C2=1, C3=2, C4=3, C5=4, C6=5, N1=6, N2=7, O1=8, O2=9, S1=10)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--species', required=True, choices=sorted(SPECIES))
    ap.add_argument('--xyz', required=True, help='reference geometry (xyz)')
    ap.add_argument('--set-angle', nargs=4, metavar=('A', 'B', 'C', 'DEG'), required=True,
                    help='set the A-B-C angle to DEG degrees (atom labels, or --index positions)')
    ap.add_argument('--index', action='store_true', help='--set-angle takes 1-based indices')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    path = pathlib.Path(a.xyz)
    path = path if path.is_absolute() else ROOT / path
    mol = from_dft_xyz(path, SPECIES[a.species][0], a.species)
    i, j, k = (int(x) - 1 for x in a.set_angle[:3]) if a.index \
        else (LABELS[x] for x in a.set_angle[:3])
    deg = float(a.set_angle[3])
    conf = mol.GetConformer()
    before = T.GetAngleDeg(conf, i, j, k)
    T.SetAngleDeg(conf, i, j, k, deg)
    after = T.GetAngleDeg(conf, i, j, k)
    out = pathlib.Path(a.out)
    out = out if out.is_absolute() else ROOT / out
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(mol.GetNumAtoms()),
             f'{a.species} from {path.name}, angle {"-".join(a.set_angle[:3])} '
             f'{before:.2f} -> {after:.2f} deg (diagnostic geometry, nothing else changed)']
    for at in mol.GetAtoms():
        p = conf.GetAtomPosition(at.GetIdx())
        lines.append(f'{at.GetSymbol():2s} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}')
    out.write_text('\n'.join(lines) + '\n')
    # report every internal coordinate that moved, so "one coordinate changed" is a measurement
    ref = from_dft_xyz(path, SPECIES[a.species][0], a.species)
    rc = ref.GetConformer()
    moved_b = max(abs(T.GetBondLength(conf, b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                      - T.GetBondLength(rc, b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
                  for b in mol.GetBonds())
    print(f'wrote {out}  angle {before:.2f} -> {after:.2f} deg, largest bond-length change '
          f'{moved_b:.4f} A')
    return 0


if __name__ == '__main__':
    sys.exit(main())
