"""Evidence for the protonation, conformer, water and start-pose decisions of the MD system (stage S3).

    /root/cqc-sfa-tools/dock/bin/python md/system_audit.py

Reads only: the deposited model (dock/pdb/6RL9.cif.gz), the docking receptor and its reports
(dock/T3_6RL9/), the docking outputs (dock/T3_6RL9/runs/).  Writes md/system/audit/:

    deposited.json        what the entry holds (hetero groups, hydrogens, termini, zinc geometry)
    altlocs.csv           every residue with alternate conformers, the label taken here and by the docking
    histidines.json       per histidine and conformer: deposited N-H, partners of ND1 and NE2, distance to Zn
    glu106.json           Glu106 / Glu117 / Thr199 geometry, carboxylate bond lengths
    propka.json           PROPKA 3.5.1 on the docking's heavy-atom protein, without and with the zinc
    propka_without_zn.pka, propka_with_zn.pka
    amides.json           Asn/Gln amide orientation check
    waters.csv            the 216 crystal waters
    poses_ad4zn.csv       every mode of the AutoDock4Zn runs of neutral sulfanilamide in 6RL9
    poses_summary.json    near-native / ring-rotated modes per run, ranges of the orientation distances

No decision is taken here; md/build_system.py holds the decisions.  Running it again gives the same
files (no random numbers; PROPKA is deterministic).
"""
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

import gemmi
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import system_common as C  # noqa: E402

AUD = C.OUT / 'audit'
POLAR = ('N', 'O', 'S', 'Zn', 'Hg')
AD4ZN_JOBS = ['sfa_neutral__ad4zn__prod', 'sfa_neutral__ad4zn', 'sfa_neutral__ad4zn__s20',
              'sfa_neutral__ad4zn__etkdg', 'sfa_neutral__ad4zn__csn102', 'sfa_neutral__ad4zn__csn109']


def r2(x, n=2):
    return round(float(x), n)


def tag(cra, sym=False):
    a = cra.atom
    return (f'{cra.residue.name}{cra.residue.seqid.num}.{a.name}'
            + (f'[{C.alt(a)}]' if C.alt(a) != '-' else '')
            + (f'(occ {a.occ:.2f})' if a.occ < 0.995 else '') + ('*' if sym else ''))


class Model:
    def __init__(self):
        self.st, self.chain = C.read_deposited()
        self.model = self.st[0]
        self.ns = gemmi.NeighborSearch(self.model, self.st.cell, 6).populate()
        self.choice, self.groups = C.choose_conformers(self.chain)

    def find(self, num, name, label=None):
        for r in self.chain:
            if r.seqid.num == num:
                for a in r:
                    if a.name == name and (C.alt(a) == '-' or C.alt(a) == (label or self.choice.get(num, 'A'))):
                        return r, a
        raise KeyError((num, name, label))

    def partners(self, res, atom, cut=3.5, elements=POLAR, with_xyz=False):
        """Atoms within cut of `atom` (same label or unlabelled), symmetry mates marked '*'."""
        out = []
        for mk in self.ns.find_atoms(atom.pos, atom.altloc, radius=cut):
            cra = mk.to_cra(self.model)
            if cra.residue is res and mk.image_idx == 0:
                continue
            if elements and cra.atom.element.name not in elements:
                continue
            p = self.st.cell.find_nearest_pbc_position(atom.pos, cra.atom.pos, mk.image_idx)
            row = (r2(atom.pos.dist(p)), tag(cra, mk.image_idx != 0))
            out.append(row + (np.array(p.tolist()),) if with_xyz else row)
        return sorted(out, key=lambda t: t[:2])

    def donor_geometry(self, res, at, which):
        """If the ring nitrogen `which` (ND1 or NE2) carried a hydrogen (in the ring plane, on the outer
        bisector, 1.01 A): H...X distance and N-H...X angle to each partner within 3.5 A."""
        n = C.xyz(at[which])
        c1 = C.xyz(at['CG' if which == 'ND1' else 'CD2'])
        c2 = C.xyz(at['CE1'])
        v = (n - c1) / np.linalg.norm(n - c1) + (n - c2) / np.linalg.norm(n - c2)
        h = n + 1.01 * v / np.linalg.norm(v)
        return [dict(partner=t, N_X=d, H_X=r2(C.dist(h, p)), N_H_X=r2(C.angle(n, h, p), 0))
                for d, t, p in self.partners(res, at[which], with_xyz=True) if not t.startswith(('ZN', 'HG'))]


def deposited(m):
    st_h, ch_h = C.read_deposited(with_h=True)
    nh = sum(1 for r in ch_h for a in r if a.element.name in ('H', 'D'))
    poly = C.polymer_residues(m.chain)
    nums = [r.seqid.num for r in poly]
    ent = next(e for e in st_h.entities if e.entity_type == gemmi.EntityType.Polymer)
    het = []
    zn0 = m.find(C.ZN_SEQ, 'ZN')[1]
    lig0 = [a for r in m.chain if r.seqid.num == C.LIG_SEQ for a in r]
    for r in C.hetero(m.chain):
        near = set()
        for a in r:
            for mk in m.ns.find_atoms(a.pos, '\0', radius=4.0):
                cra = mk.to_cra(m.model)
                if cra.residue.het_flag == 'A':
                    near.add(f'{cra.residue.name}{cra.residue.seqid.num}' + ('*' if mk.image_idx else ''))
        het.append(dict(name=r.name, seq=r.seqid.num, n_atoms=len(r), labels=C.labels_of(r) or ['-'],
                        occupancy=sorted({r2(a.occ) for a in r}),
                        B_range=[r2(min(a.b_iso for a in r), 1), r2(max(a.b_iso for a in r), 1)],
                        min_dist_to_Zn_A=r2(min(a.pos.dist(zn0.pos) for a in r)),
                        min_dist_to_SAN304_A=r2(min(a.pos.dist(b.pos) for a in r for b in lig0)),
                        protein_residues_within_4A=sorted(near, key=lambda t: (t.endswith('*'), int(''.join(c for c in t if c.isdigit())))),
                        note='* = residue of a neighbouring molecule of the crystal'))
    lig = [(a.name, C.alt(a), r2(a.occ), r2(a.b_iso, 1)) for r in m.chain if r.seqid.num == C.LIG_SEQ for a in r]
    zn = m.find(C.ZN_SEQ, 'ZN')[1]
    site = [(94, 'NE2'), (96, 'NE2'), (119, 'ND1'), (C.LIG_SEQ, 'N2')]
    zd = {f'{n}.{nm}': r2(zn.pos.dist(m.find(n, nm)[1].pos), 3) for n, nm in site}
    za = {}
    for i in range(4):
        for j in range(i + 1, 4):
            a, b = m.find(*site[i])[1], m.find(*site[j])[1]
            za[f'{site[i][0]}.{site[i][1]}-Zn-{site[j][0]}.{site[j][1]}'] = r2(
                C.angle(C.xyz(a), C.xyz(zn), C.xyz(b)), 1)
    first, last = ch_h[0], [r for r in ch_h if r.seqid.num == nums[-1]][0]
    w = C.waters(m.chain)
    return dict(
        entry='6RL9', file=str(C.CIF.relative_to(C.REPO)), sha256=C.sha256(C.CIF),
        resolution_A=st_h.resolution, spacegroup=st_h.spacegroup_hm,
        deposited_sequence_length=len(ent.full_sequence),
        deposited_sequence_start=' '.join(ent.full_sequence[:9]),
        modelled_residues=dict(first=f'{poly[0].name}{nums[0]}', last=f'{poly[-1].name}{nums[-1]}', n=len(poly)),
        numbering_gaps=[(a, b) for a, b in zip(nums, nums[1:]) if b - a != 1],
        peptide_bond_across_gap_A={f'{a}-{b}': r2(m.find(a, 'C')[1].pos.dist(m.find(b, 'N')[1].pos), 3)
                                   for a, b in zip(nums, nums[1:]) if b - a != 1},
        hydrogens_deposited=nh,
        hydrogens_on_first_residue_N=[a.name for a in first if a.element.name == 'H' and a.name in ('H', 'H1', 'H2', 'H3')],
        last_residue_heavy_atoms=sorted({a.name for a in last if a.element.name != 'H'}),
        hydrogens_on_ligand_304=sum(1 for r in ch_h if r.seqid.num == C.LIG_SEQ for a in r if a.element.name == 'H'),
        hydrogens_on_waters=sum(1 for r in ch_h if r.name == 'HOH' for a in r if a.element.name == 'H'),
        hetero_groups=het, ligand_304_atoms=lig,
        n_waters=len(w), n_water_sites_with_label=sum(1 for r in w if C.labels_of(r)),
        zinc=dict(B=r2(zn.b_iso, 1), distances_A=zd, angles_deg=za,
                  Zn_O1=r2(zn.pos.dist(m.find(C.LIG_SEQ, 'O1')[1].pos), 3),
                  Zn_O2=r2(zn.pos.dist(m.find(C.LIG_SEQ, 'O2')[1].pos), 3)),
    )


def altlocs(m):
    dock = C.docking_conformers()
    lig = [a for r in m.chain if r.seqid.num == C.LIG_SEQ for a in r]
    zn = m.find(C.ZN_SEQ, 'ZN')[1]
    rows = []
    for r in C.polymer_residues(m.chain):
        labs = C.labels_of(r)
        if not labs:
            continue
        by = {L: {a.name: a for a in r if C.alt(a) == L} for L in labs}
        disp = []
        for L in labs[1:]:
            common = set(by[labs[0]]) & set(by[L])
            bb = max((by[labs[0]][n].pos.dist(by[L][n].pos) for n in C.BACKBONE if n in common), default=0)
            sc = max((by[labs[0]][n].pos.dist(by[L][n].pos) for n in common if n not in C.BACKBONE), default=0)
            disp.append(f'{labs[0]}/{L}: backbone {bb:.2f}, side chain {sc:.2f}')
        grp = next(g for g in m.groups if f'{r.name}{r.seqid.num}' in g['residues'])
        rows.append(dict(
            residue=f'{r.name}{r.seqid.num}', labels='/'.join(labs),
            occupancy='/'.join(f'{next(iter(by[L].values())).occ:.2f}' for L in labs),
            heavy_atoms='/'.join(str(len(by[L])) for L in labs),
            expected_heavy_atoms=len(C.HEAVY[r.name]),
            max_displacement_A='; '.join(disp),
            group='+'.join(grp['residues']), decided_by=grp['decided_by'],
            coupling=' | '.join(grp['coupling']),
            taken_md=m.choice[r.seqid.num], taken_docking=dock.get(r.seqid.num, ''),
            differs_from_docking=m.choice[r.seqid.num] != dock.get(r.seqid.num, ''),
            min_dist_to_SAN304_A=r2(min(a.pos.dist(b.pos) for a in r for b in lig)),
            min_dist_to_Zn_A=r2(min(a.pos.dist(zn.pos) for a in r)),
        ))
    return rows


def histidines(m):
    st_h, ch_h = C.read_deposited(with_h=True)
    zn = m.find(C.ZN_SEQ, 'ZN')[1]
    out = []
    for r in C.polymer_residues(m.chain):
        if r.name != 'HIS':
            continue
        rh = [q for q in ch_h if q.seqid.num == r.seqid.num][0]
        for L in (C.labels_of(r) or ['-']):
            at = {a.name: a for a in r if C.alt(a) in ('-', L)}
            dep = sorted(a.name for a in rh if a.name in ('HD1', 'HE2') and C.alt(a) in ('-', L))
            chi1 = C.dihedral(*(C.xyz(at[n]) for n in ('N', 'CA', 'CB', 'CG')))
            chi2 = C.dihedral(*(C.xyz(at[n]) for n in ('CA', 'CB', 'CG', 'ND1')))
            out.append(dict(
                residue=f'HIS{r.seqid.num}', label=L, occupancy=r2(at['CA'].occ),
                taken=(L == '-' or m.choice.get(r.seqid.num) == L),
                deposited_NH=dep, B_ring=r2(np.mean([at[n].b_iso for n in ('ND1', 'CD2', 'CE1', 'NE2')]), 1),
                chi1=r2(chi1, 1), chi2=r2(chi2, 1),
                Zn_ND1=r2(zn.pos.dist(at['ND1'].pos)), Zn_NE2=r2(zn.pos.dist(at['NE2'].pos)),
                ND1_partners=m.partners(r, at['ND1']), NE2_partners=m.partners(r, at['NE2']),
                if_H_on_ND1=m.donor_geometry(r, at, 'ND1'), if_H_on_NE2=m.donor_geometry(r, at, 'NE2'),
            ))
    return out


def glu106(m):
    st_h, ch_h = C.read_deposited(with_h=True)
    P = lambda n, nm: C.xyz(m.find(n, nm)[1])  # noqa: E731
    hg1 = C.xyz([a for r in ch_h if r.seqid.num == 199 for a in r if a.name == 'HG1'][0])
    zn = P(C.ZN_SEQ, 'ZN')
    out = dict(
        note='distances in A, angles in degrees; deposited model, hydrogens are the riding hydrogens of the entry',
        Glu106=dict(OE1_partners=m.partners(*m.find(106, 'OE1')), OE2_partners=m.partners(*m.find(106, 'OE2')),
                    CD_OE1=r2(C.dist(P(106, 'CD'), P(106, 'OE1')), 3), CD_OE2=r2(C.dist(P(106, 'CD'), P(106, 'OE2')), 3),
                    Zn_OE1=r2(C.dist(zn, P(106, 'OE1'))), Zn_OE2=r2(C.dist(zn, P(106, 'OE2'))),
                    B=[r2(m.find(106, n)[1].b_iso, 1) for n in ('OE1', 'OE2')],
                    deposited_carboxyl_H=[a.name for r in ch_h if r.seqid.num == 106 for a in r if a.name in ('HE1', 'HE2')]),
        Glu117=dict(OE1_partners=m.partners(*m.find(117, 'OE1')), OE2_partners=m.partners(*m.find(117, 'OE2')),
                    CD_OE1=r2(C.dist(P(117, 'CD'), P(117, 'OE1')), 3), CD_OE2=r2(C.dist(P(117, 'CD'), P(117, 'OE2')), 3),
                    Zn_OE1=r2(C.dist(zn, P(117, 'OE1'))), Zn_OE2=r2(C.dist(zn, P(117, 'OE2')))),
        Thr199=dict(OG1_partners=m.partners(*m.find(199, 'OG1')), N_partners=m.partners(*m.find(199, 'N'), cut=3.3),
                    OG1_OE1=r2(C.dist(P(199, 'OG1'), P(106, 'OE1'))), OG1_OE2=r2(C.dist(P(199, 'OG1'), P(106, 'OE2'))),
                    OG1_N2=r2(C.dist(P(199, 'OG1'), P(C.LIG_SEQ, 'N2'))),
                    deposited_HG1=dict(HG1_OE1=r2(C.dist(hg1, P(106, 'OE1'))),
                                       OG1_HG1_OE1=r2(C.angle(P(199, 'OG1'), hg1, P(106, 'OE1')), 1),
                                       HG1_N2=r2(C.dist(hg1, P(C.LIG_SEQ, 'N2'))),
                                       OG1_HG1_N2=r2(C.angle(P(199, 'OG1'), hg1, P(C.LIG_SEQ, 'N2')), 1)),
                    N_O1=r2(C.dist(P(199, 'N'), P(C.LIG_SEQ, 'O1'))), N_O2=r2(C.dist(P(199, 'N'), P(C.LIG_SEQ, 'O2')))),
        water415=dict(partners=m.partners(*m.find(415, 'O'))),
    )
    # carboxylate bond lengths of every fully occupied Asp/Glu: does the refinement let them differ?
    diffs = []
    for r in C.polymer_residues(m.chain):
        if r.name not in ('ASP', 'GLU') or C.labels_of(r):
            continue
        c, o1, o2 = ('CG', 'OD1', 'OD2') if r.name == 'ASP' else ('CD', 'OE1', 'OE2')
        at = {a.name: a for a in r}
        if not all(n in at for n in (c, o1, o2)):
            continue
        d1, d2 = at[c].pos.dist(at[o1].pos), at[c].pos.dist(at[o2].pos)
        diffs.append((abs(d1 - d2), f'{r.name}{r.seqid.num}', r2(d1, 3), r2(d2, 3)))
    diffs.sort(reverse=True)
    v = np.array([d[0] for d in diffs])
    out['carboxylate_CO_difference'] = dict(
        n=len(diffs), mean_abs_diff=r2(v.mean(), 3), median=r2(np.median(v), 3), max=r2(v.max(), 3),
        largest=[dict(residue=d[1], d1=d[2], d2=d[3]) for d in diffs[:5]],
        Glu106_abs_diff=r2(abs(out['Glu106']['CD_OE1'] - out['Glu106']['CD_OE2']), 3),
        rank_of_Glu106=1 + sum(1 for d in diffs if d[0] > abs(out['Glu106']['CD_OE1'] - out['Glu106']['CD_OE2']) + 1e-9),
        note='a protonated carboxyl group has C=O near 1.21 and C-OH near 1.31 A; geometry restraints of the refinement pull both towards 1.25 A, so this is weak evidence')
    return out


def propka():
    """PROPKA on the heavy-atom protein the docking used, without the zinc (what the docking did) and
    with the zinc as a +2 ion (PROPKA knows ZN)."""
    heavy = (C.DOCK / 'T3_6RL9' / 'fixed_heavy.pdb').read_text().splitlines()
    heavy = [ln for ln in heavy if ln.startswith('ATOM')]
    znl = [ln for ln in C.DOCK_RECEPTOR.read_text().splitlines() if ln.startswith('HETATM') and ' ZN ' in ln]
    res = {}
    with tempfile.TemporaryDirectory() as tmp:
        for name, lines in (('without_zn', heavy), ('with_zn', heavy + znl)):
            f = os.path.join(tmp, f'{name}.pdb')
            open(f, 'w').write('\n'.join(lines + ['END']) + '\n')
            p = subprocess.run([sys.executable, '-m', 'propka', f'{name}.pdb'], cwd=tmp, capture_output=True, text=True)
            if p.returncode != 0:
                raise SystemExit('PROPKA failed: ' + p.stderr[-500:])
            txt = open(os.path.join(tmp, f'{name}.pka')).read()
            shutil.copy(os.path.join(tmp, f'{name}.pka'), AUD / f'propka_{name}.pka')
            block = txt.split('SUMMARY OF THIS PREDICTION')[1].split('-----')[0]
            pk = {}
            for ln in block.splitlines():
                mm = re.match(r'\s+([A-Z+\-]+)\s+(\d+)\s+A\s+(-?[\d.]+)\s+(-?[\d.]+)', ln)
                if mm:
                    pk[f'{mm.group(1)}{mm.group(2)}'] = (float(mm.group(3)), float(mm.group(4)))
            res[name] = pk
    import propka as _p
    tab, flagged = {}, []
    for k in res['without_zn']:
        a, model = res['without_zn'][k]
        b = res['with_zn'].get(k, (None, None))[0]
        tab[k] = dict(without_zn=a, with_zn=b, model_pKa=model)
        kind = re.match(r'[A-Z+\-]+', k).group(0)
        for label, v in (('without_zn', a), ('with_zn', b)):
            if v is None:
                continue
            unusual = (v > 6.4) if kind in ('ASP', 'GLU', 'C-') else (v < 8.4) if kind in ('TYR', 'CYS', 'LYS', 'ARG', 'N+') \
                else (v > 6.4) if kind == 'HIS' else False
            if unusual:
                flagged.append(dict(group=k, run=label, pKa=v))
    return dict(propka_version=getattr(_p, '__version__', '3.5.1'), pH=7.4,
                input='dock/T3_6RL9/fixed_heavy.pdb (the heavy-atom protein of the docking), zinc from dock/T3_6RL9/receptor.pdb',
                flag_rule='acids with pKa > 6.4, bases (Lys, Arg, Tyr, Cys, N-terminus) with pKa < 8.4, His with pKa > 6.4 (within one unit of pH 7.4 or on the unusual side)',
                flagged=flagged, pKa=tab)


def amides(m):
    """Asn/Gln: is the deposited orientation of the amide consistent with its neighbours?  The oxygen
    should not face a carboxylate or carbonyl oxygen only, the nitrogen should not face an N-H donor only."""
    lig = [a for r in m.chain if r.seqid.num == C.LIG_SEQ for a in r]
    acc_only = {'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OXT'}
    out = []
    for r in C.polymer_residues(m.chain):
        if r.name not in ('ASN', 'GLN'):
            continue
        L = m.choice.get(r.seqid.num, '-')
        at = {a.name: a for a in r if C.alt(a) in ('-', L)}
        o, n = ('OD1', 'ND2') if r.name == 'ASN' else ('OE1', 'NE2')
        if o not in at or n not in at:
            continue
        po, pn = m.partners(r, at[o], cut=3.3, elements=('N', 'O')), m.partners(r, at[n], cut=3.3, elements=('N', 'O'))

        def kinds(pl):
            k = []
            for d, t in pl:
                if t.endswith('*'):
                    continue
                mm = re.match(r'([A-Z0-9]{2,3}?)(\d+)\.([A-Z0-9]+)', t)
                rn, an = mm.group(1), mm.group(3)
                if rn == 'HOH':
                    k.append('water')
                elif an in acc_only and not (rn in ('ASN', 'GLN') and an in ('OD1', 'OE1')):
                    k.append('acceptor')
                elif an in ('N', 'NZ', 'NE', 'NH1', 'NH2', 'NE1'):
                    k.append('donor')
                else:
                    k.append('either')
            return k
        ko, kn = kinds(po), kinds(pn)
        bad = ('acceptor' in ko and 'donor' not in ko and 'water' not in ko and 'either' not in ko) or \
              ('donor' in kn and 'acceptor' not in kn and 'water' not in kn and 'either' not in kn)
        out.append(dict(residue=f'{r.name}{r.seqid.num}', label=L, O_atom=o, N_atom=n,
                        O_partners=po, N_partners=pn, B_O=r2(at[o].b_iso, 1), B_N=r2(at[n].b_iso, 1),
                        min_dist_to_SAN304_A=r2(min(a.pos.dist(b.pos) for a in at.values() for b in lig)),
                        suspicious=bool(bad)))
    return out


def water_table(m):
    zn = m.find(C.ZN_SEQ, 'ZN')[1]
    lig = [a for r in m.chain if r.seqid.num == C.LIG_SEQ for a in r]
    kept, kept_tag = [], []
    for r in C.polymer_residues(m.chain):
        for a in C.conformer_atoms(r, m.choice.get(r.seqid.num, '-')):
            kept.append(a)
            kept_tag.append(f'{r.name}{r.seqid.num}.{a.name}')
    kept_pos = np.array([C.xyz(a) for a in kept])
    rejected = [(r, a) for r in C.polymer_residues(m.chain) for a in r
                if C.alt(a) not in ('-', m.choice.get(r.seqid.num, '-'))]
    het = [(r, a) for r in C.hetero(m.chain) if r.seqid.num not in (C.ZN_SEQ, C.LIG_SEQ) for a in r]
    rows = []
    for r in C.waters(m.chain):
        for a in r:
            p = C.xyz(a)
            d = np.linalg.norm(kept_pos - p, axis=1)
            i = int(d.argmin())
            dr = min(((a.pos.dist(b.pos), f'{q.name}{q.seqid.num}.{b.name}[{C.alt(b)}]') for q, b in rejected), default=(99, ''))
            dh = min(((a.pos.dist(b.pos), f'{q.name}{q.seqid.num}.{b.name}') for q, b in het), default=(99, ''))
            dl = min((a.pos.dist(b.pos), b.name) for b in lig)
            rows.append(dict(
                water=r.seqid.num, label=C.alt(a), occupancy=r2(a.occ), B=r2(a.b_iso, 1),
                x=r2(p[0], 3), y=r2(p[1], 3), z=r2(p[2], 3),
                min_dist_protein_taken_A=r2(d[i]), nearest_protein_atom=kept_tag[i],
                min_dist_rejected_conformer_A=r2(dr[0]), nearest_rejected=dr[1] if dr[0] < 4 else '',
                min_dist_removed_hetero_A=r2(dh[0]), nearest_removed_hetero=dh[1] if dh[0] < 4 else '',
                dist_Zn_A=r2(a.pos.dist(zn.pos)), min_dist_SAN304_A=r2(dl[0]), nearest_SAN304_atom=dl[1],
                polar_partners=' '.join(f'{t}:{d_}' for d_, t in m.partners(r, a, cut=3.3, elements=('N', 'O'))),
            ))
    return rows


def pose_tables(m):
    site = C.site_atoms(m.chain, m.choice)
    crystal = {n: C.xyz(m.find(C.LIG_SEQ, n)[1]) for n in C.SAN_HEAVY}
    rows = []
    for job in AD4ZN_JOBS:
        d = C.DOCK / 'T3_6RL9' / 'runs' / job
        for f in sorted(d.glob('seed*/out.pdbqt')):
            seed = int(f.parent.name[4:])
            for i, mode in enumerate(C.read_pdbqt_modes(f), 1):
                rmsd, M, mp = C.name_docked_atoms(mode, crystal)
                mt = C.pose_metrics(M, site)
                rows.append(dict(job=job, seed=seed, mode=i, score=mode['score'], rmsd_A=round(rmsd, 3),
                                 **{k: round(v, 3) for k, v in mt.items()}))
    summ = {}
    for job in AD4ZN_JOBS:
        g = [r for r in rows if r['job'] == job]
        ok = [r for r in g if r['Zn_N2'] <= 2.3 and r['N2_OG1'] <= 3.3]
        near = [r for r in ok if r['rmsd_A'] <= 1.0]
        rot = [r for r in ok if r['rmsd_A'] > 2.0]
        # the ring-rotated orientation of the paper: aniline N at Asn67 OD1 (within 3.5 A)
        rot67 = [r for r in rot if r['N1_Asn67OD1'] <= 3.5]
        top = [r for r in g if r['mode'] == 1]

        def rng(rs, k):
            return [min(r[k] for r in rs), max(r[k] for r in rs)] if rs else None
        gaps, ranks = [], []
        for seed in sorted({r['seed'] for r in g}):
            nn = [r for r in near if r['seed'] == seed]
            rr = [r for r in rot67 if r['seed'] == seed]
            if nn and rr:
                gaps.append(round(min(r['score'] for r in rr) - min(r['score'] for r in nn), 3))
                ranks.append(min(r['mode'] for r in rr))
        summ[job] = dict(
            seeds=len({r['seed'] for r in g}), modes=len(g), modes_passing_site_check=len(ok),
            near_native_modes=len(near), modes_beyond_2A=len(rot), rotated_modes_at_Asn67=len(rot67),
            top_poses_near_native=sum(1 for r in top if r['rmsd_A'] <= 1.0),
            top_poses_rotated=sum(1 for r in top if r['rmsd_A'] > 2.0),
            near_native=dict(rmsd=rng(near, 'rmsd_A'), N1_Phe131CZ=rng(near, 'N1_Phe131CZ'), N1_Asn67OD1=rng(near, 'N1_Asn67OD1')),
            rotated_at_Asn67=dict(rmsd=rng(rot67, 'rmsd_A'), score=rng(rot67, 'score'),
                                  N1_Phe131CZ=rng(rot67, 'N1_Phe131CZ'), N1_Asn67OD1=rng(rot67, 'N1_Asn67OD1')),
            seeds_with_both_orientations=len(gaps),
            score_best_rotated_minus_best_near_native_per_seed=[min(gaps), max(gaps)] if gaps else None,
            rank_of_first_rotated_mode_per_seed=sorted(set(ranks)),
            rotated_top_poses=dict(rmsd=rng([r for r in top if r['rmsd_A'] > 2.0], 'rmsd_A'),
                                   N1_Phe131CZ=rng([r for r in top if r['rmsd_A'] > 2.0], 'N1_Phe131CZ'),
                                   N1_Asn67OD1=rng([r for r in top if r['rmsd_A'] > 2.0], 'N1_Asn67OD1')),
        )
    crystal_m = {k: round(v, 3) for k, v in C.pose_metrics(crystal, site).items()}
    return rows, dict(
        note='site check as in the paper: Zn-N(sulfonamide) <= 2.3 A and Thr199 OG1...N <= 3.3 A; near-native: heavy-atom RMSD <= 1.0 A; rotated: > 2.0 A; distances to the deposited coordinates of 6RL9',
        crystal_pose=crystal_m, jobs=summ)


def write_csv(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main():
    AUD.mkdir(parents=True, exist_ok=True)
    m = Model()
    dump = lambda name, obj: (AUD / name).write_text(json.dumps(obj, indent=1) + '\n')  # noqa: E731
    dump('deposited.json', deposited(m))
    write_csv(AUD / 'altlocs.csv', altlocs(m))
    dump('histidines.json', histidines(m))
    dump('glu106.json', glu106(m))
    dump('propka.json', propka())
    dump('amides.json', amides(m))
    write_csv(AUD / 'waters.csv', water_table(m))
    rows, summ = pose_tables(m)
    write_csv(AUD / 'poses_ad4zn.csv', rows)
    dump('poses_summary.json', summ)
    print('audit written to', AUD.relative_to(C.REPO))
    return 0


if __name__ == '__main__':
    sys.exit(main())
