"""Flip one neutral histidine's tautomer in an already prepared receptor, changing nothing else.

    python dock/his_tautomer_variant.py --target T3_6RL9 --resnum 64 --to HID
    # writes dock/T3_6RL9_his64HID/{receptor.pdbqt,receptor.pdb,box.json,ref_*.sdf,variant.json}

Why not simply re-run prep_receptor.py with a different tautomer.  PDBFixer's rebuild and OpenMM's
hydrogen placement are not bit-reproducible, and the restrained minimisation moves free atoms by up to
several angstrom, so a freshly prepared receptor differs from the one the reported docking used in many
places at once.  A sensitivity test has to change one thing.  This script therefore edits the existing
receptor: it moves the single imidazole proton between NE2 and ND1 and swaps the two nitrogens'
AutoDock types, which is exactly what the HIE/HID distinction is.

    HIE: proton on NE2  ->  NE2 = N (donor, has HD), ND1 = NA (acceptor)
    HID: proton on ND1  ->  ND1 = N (donor, has HD), NE2 = NA (acceptor)

The proton is placed in the ring plane along the external bisector of the C-N-C angle at 1.01 A, which
is where Modeller would put it.  Vina's scoring function has no electrostatics, so only the types and
the proton position matter to it; the partial charges are swapped with the nitrogens anyway so that the
file stays self-consistent for AutoDock4 as well.  Nothing else in the file is touched, and
variant.json records the before/after types, the proton geometry and the sha256 of both receptors.

Scope: neutral His only (HIE <-> HID).  It refuses a residue that has both HD1 and HE2 (that is HIP,
a charged histidine, and flipping it is not a tautomer change).
"""
import argparse
import hashlib
import json
import pathlib
import shutil
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

DOCK = pathlib.Path(__file__).resolve().parent
NH = 1.01          # A, N-H bond length used by the Amber templates


def parse(path):
    return path.read_text().splitlines()


def coords(line):
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def set_coords(line, xyz):
    return f'{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}'


def rename_atom(line, name):
    """Write an atom name into columns 13-16 the way the PDB format wants it: a name of three
    characters or fewer starts in column 14, so that column 13-14 still hold the element."""
    field = f'{name:<4s}' if len(name) >= 4 else f' {name:<3s}'
    return f'{line[:12]}{field}{line[16:]}'


def place_h(n, c1, c2):
    """In-plane external bisector of C-N-C at NH from the nitrogen."""
    v = (n - c1) / np.linalg.norm(n - c1) + (n - c2) / np.linalg.norm(n - c2)
    return n + NH * v / np.linalg.norm(v)


def residue_atoms(lines, resnum, chain, resname='HIS'):
    idx = {}
    for i, ln in enumerate(lines):
        if not ln.startswith(('ATOM', 'HETATM')):
            continue
        if ln[17:20].strip() == resname and int(ln[22:26]) == resnum and ln[21] == chain:
            idx[ln[12:16].strip()] = i
    return idx


def flip(lines, resnum, chain, to, pdbqt):
    """Return (new lines, record).  `pdbqt` switches on the AutoDock type/charge columns."""
    idx = residue_atoms(lines, resnum, chain)
    if not idx:
        raise SystemExit(f'His{resnum} chain {chain} not found')
    have_hd1, have_he2 = 'HD1' in idx, 'HE2' in idx
    if have_hd1 and have_he2:
        raise SystemExit(f'His{resnum} carries both HD1 and HE2 (HIP): not a tautomer flip')
    frm = 'HID' if have_hd1 else 'HIE'
    if frm == to:
        raise SystemExit(f'His{resnum} is already {to}')
    hname_old, hname_new = ('HD1', 'HE2') if frm == 'HID' else ('HE2', 'HD1')
    n_old, n_new = ('ND1', 'NE2') if frm == 'HID' else ('NE2', 'ND1')
    for nm in (n_old, n_new, 'CG', 'CD2', 'CE1', hname_old):
        if nm not in idx:
            raise SystemExit(f'His{resnum}: atom {nm} missing, cannot flip')
    out = list(lines)
    # the proton moves to the other nitrogen, in the ring plane
    nn = coords(lines[idx[n_new]])
    c1 = coords(lines[idx['CG' if n_new == 'ND1' else 'CD2']])
    c2 = coords(lines[idx['CE1']])
    hpos = place_h(nn, c1, c2)
    h_line = set_coords(rename_atom(lines[idx[hname_old]], hname_new), hpos)
    out[idx[hname_old]] = h_line
    rec = dict(residue=f'{chain}:HIS{resnum}', frm=frm, to=to,
               proton_renamed=[hname_old, hname_new],
               proton_moved_A=round(float(np.linalg.norm(hpos - coords(lines[idx[hname_old]]))), 3),
               new_N_H_distance_A=round(float(np.linalg.norm(hpos - nn)), 3))
    if pdbqt:
        # swap the two nitrogens' AutoDock types (N = donor, NA = acceptor) and their charges
        def typ(i):
            return out[i][77:79].strip()

        def retype(i, t, q):
            return f'{out[i][:70]}{q:6.3f} {t:<2s}'
        q_old, q_new = float(out[idx[n_old]][70:76]), float(out[idx[n_new]][70:76])
        t_old, t_new = typ(idx[n_old]), typ(idx[n_new])
        rec['types_before'] = {n_old: t_old, n_new: t_new}
        out[idx[n_old]] = retype(idx[n_old], t_new, q_new)
        out[idx[n_new]] = retype(idx[n_new], t_old, q_old)
        rec['types_after'] = {n_old: t_new, n_new: t_old}
        rec['charges_swapped'] = {n_old: q_new, n_new: q_old}
    return out, rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', required=True, help='an existing prepared receptor directory name')
    ap.add_argument('--resnum', type=int, required=True)
    ap.add_argument('--chain', default=None, help='default: the chain of the target definition')
    ap.add_argument('--to', required=True, choices=['HID', 'HIE'])
    a = ap.parse_args()
    from targets import TARGETS
    src = DOCK / a.target
    chain = a.chain or TARGETS[a.target.split('_seed')[0]]['chain']
    dst = DOCK / f'{a.target}_his{a.resnum}{a.to}'
    dst.mkdir(parents=True, exist_ok=True)
    record = dict(source=a.target, variant=dst.name, chain=chain, resnum=a.resnum, to=a.to,
                  note='single-residue tautomer flip of an existing receptor; nothing else changed')
    for name, is_pdbqt in (('receptor.pdbqt', True), ('receptor.pdb', False)):
        lines = parse(src / name)
        new, rec = flip(lines, a.resnum, chain, a.to, is_pdbqt)
        (dst / name).write_text('\n'.join(new) + '\n')
        record[name] = rec
        record[f'{name}_sha256'] = dict(
            source=hashlib.sha256((src / name).read_bytes()).hexdigest(),
            variant=hashlib.sha256((dst / name).read_bytes()).hexdigest())
        assert len(new) == len(lines)
    # The box, the Vina config and the target's own co-crystal ligand come along; ref_cross_*.sdf are
    # analyse.py's own cross-docking products and are left behind so it regenerates them in this frame.
    own_refs = [p.name for p in src.glob('ref_*.sdf') if not p.name.startswith('ref_cross_')]
    for extra in ['box.json', 'vina_config.txt', *own_refs]:
        if (src / extra).exists():
            shutil.copy2(src / extra, dst / extra)
    (dst / 'variant.json').write_text(json.dumps(record, indent=2) + '\n')
    r = record['receptor.pdbqt']
    print(f"{dst.name}: {r['residue']} {r['frm']} -> {r['to']}, proton {r['proton_renamed'][0]} -> "
          f"{r['proton_renamed'][1]} moved {r['proton_moved_A']} A (N-H {r['new_N_H_distance_A']} A), "
          f"types {r['types_before']} -> {r['types_after']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
