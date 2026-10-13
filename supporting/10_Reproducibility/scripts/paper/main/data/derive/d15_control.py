"""d15: the control receptor (E. coli DNA gyrase B, PDB 4DUH) re-derived from the raw docking outputs, without
dock/control_gyrb.py, analyse.py or production_analysis.py (derive/ convention: own parsers, own geometry).

Scores   the first 'REMARK VINA RESULT' of every seedNN/out.pdbqt (the top pose); mean and SD (ddof 1) per job.
RMSD     heavy atoms of the top pose in the receptor frame (no superposition) against the crystal ligand
         (dock/C1_4DUH/ref_RLI.sdf); symmetry by enumerating the element-preserving graph isomorphisms on
         distance-perceived graphs, as in d09.
Site     a pose is in the ATP site when its heavy-atom centroid lies within 4.0 A of a heavy atom of the crystal
         inhibitor; the polar contact to Asp73 is the shortest distance from a ligand N or O to OD1/OD2.
Clorobiocin  its 1KZN pose is carried into the 4DUH frame here by a route of our own (the two sequences aligned with
         difflib, Kabsch fit of the paired CA atoms of the deposited files) and compared with the reference file
         the docking analysis wrote; its overlap with the prepared receptor is measured on that file, pair by pair,
         with the residues of 4DUH that have no counterpart in the 1KZN model marked.  Clorobiocin is not docked.
Size     Pearson r and least-squares slope of the mean score against the number of heavy atoms, over the species
         docked into 4DUH; and r over the eight species of the production matrix in every receptor.
Targets  the same species in 3TZF, 7S2K, 6RL9 and 3HS4: the raw-output scores of d09.
Loop     distance of the residues PDBFixer built (109-116) from the crystal inhibitor and from every top pose.
"""
import difflib
import glob
import gzip
import statistics as stt

import gemmi
from rdkit import Chem
from rdkit.Chem import rdchem
from common import *

C1, C1B, C1W, C1N = 'C1_4DUH', 'C1_4DUH_seed12345', 'C1_4DUH_wat', 'C1_4DUH_noloop'
SITE_CUT = 4.0
AD2EL = {'A': 'C', 'C': 'C', 'N': 'N', 'NA': 'N', 'NS': 'N', 'OA': 'O', 'OS': 'O', 'S': 'S', 'SA': 'S', 'HD': 'H',
         'H': 'H', 'F': 'F', 'Cl': 'Cl', 'CL': 'Cl', 'P': 'P'}
RC = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'P': 1.07, 'S': 1.05, 'F': 0.57, 'Cl': 1.02}
SPECIES8 = ['sfa_neutral', 'sfa_cation', 'sfa_anion', 'smx_anion', 'smx_neutral', 'paba_anion', 'azm_neutral',
            'azm_anion']


def models(pdbqt):
    """[(score, [(el, xyz)])] for every MODEL of a Vina output, heavy atoms only."""
    out, score, atoms = [], None, []
    for line in Path(pdbqt).read_text().splitlines():
        if line.startswith('MODEL'):
            score, atoms = None, []
        elif line.startswith('REMARK VINA RESULT'):
            score = float(line.split()[3])
        elif line.startswith(('ATOM', 'HETATM')):
            t = line[77:79].strip()
            el = AD2EL.get(t, t[0])
            if el != 'H':
                atoms.append((el, np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
        elif line.startswith('ENDMDL'):
            out.append((score, atoms))
    return out


def graph(atoms):
    m = Chem.RWMol()
    for el, _ in atoms:
        a = Chem.Atom(el)
        a.SetNoImplicit(True)
        m.AddAtom(a)
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if dist(atoms[i][1], atoms[j][1]) < 1.2 * (RC[atoms[i][0]] + RC[atoms[j][0]]):
                m.AddBond(i, j, rdchem.BondType.SINGLE)
    m = m.GetMol()
    m.UpdatePropertyCache(strict=False)
    return m


def sdf_atoms(path):
    mol = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=False)
    conf = mol.GetConformer()
    return [(a.GetSymbol(), np.array(conf.GetAtomPosition(a.GetIdx()))) for a in mol.GetAtoms() if a.GetSymbol() != 'H']


def rmsd_sym(pose, ref):
    """min RMSD over all element-preserving mappings of the reference graph onto the pose graph."""
    gp, gr = graph(pose), graph(ref)
    matches = gp.GetSubstructMatches(gr, uniquify=False, maxMatches=200000)
    assert matches, 'pose and reference are not the same heavy-atom graph'
    P = np.array([p[1] for p in pose])
    R = np.array([r_[1] for r_ in ref])
    return float(min(np.sqrt(np.mean(np.sum((P[list(mt)] - R) ** 2, axis=1))) for mt in matches)), len(matches)


def receptor(target):
    rows = []
    for ln in (DOCK / target / 'receptor.pdb').read_text().splitlines():
        if ln.startswith(('ATOM', 'HETATM')) and (ln[76:78].strip() or ln[12:16].strip()[0]) != 'H':
            rows.append((int(ln[22:26]), ln[17:20].strip(), ln[12:16].strip(),
                         np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
    return rows


def run(target, tag):
    seeds = sorted(glob.glob(str(DOCK / target / 'runs' / tag / 'seed*' / 'out.pdbqt')))
    return [models(s) for s in seeds]


def score_stats(per):
    sc = [m[0][0] for m in per]
    return dict(n=len(sc), mean=stt.mean(sc), sd=stt.stdev(sc) if len(sc) > 1 else 0.0, best=min(sc), per_seed=sc)


REF = sdf_atoms(DOCK / C1 / 'ref_RLI.sdf')
REFX = np.array([a[1] for a in REF])
out = {'scores': {}, 'site': {}, 'rmsd': {}, 'whole': {}}

# ---------------------------------------------------------------------------------------- every focused job
JOBS = {}
for target in (C1, C1B, C1W, C1N):
    if not (DOCK / target / 'runs').is_dir():
        continue
    for d in sorted((DOCK / target / 'runs').glob('*')):
        if d.name.endswith('__whole') or not (d / 'results.json').exists():
            continue                                  # whole-fragment runs are read below; unfinished runs are skipped
        per = run(target, d.name)
        key = f'{target}/runs/{d.name}'
        JOBS[key] = per
        s = score_stats(per)
        s['heavy_atoms'] = len(per[0][0][1])
        out['scores'][key] = s
        rec = receptor(target)
        asp = np.array([r_[3] for r_ in rec if r_[0] == 73 and r_[2] in ('OD1', 'OD2')])
        loop = np.array([r_[3] for r_ in rec if 109 <= r_[0] <= 116])
        ins, dasp, dloop, dcen = [], [], [], []
        for m in per:
            X = np.array([a[1] for a in m[0][1]])
            polar = np.array([a[1] for a in m[0][1] if a[0] in 'NO'])
            dcen.append(float(np.linalg.norm(REFX - X.mean(0), axis=1).min()))
            ins.append(dcen[-1] <= SITE_CUT)
            dasp.append(float(np.linalg.norm(asp[:, None] - polar[None], axis=2).min()))
            if len(loop):
                dloop.append(float(np.linalg.norm(loop[:, None] - X[None], axis=2).min()))
        out['site'][key] = dict(in_site_fraction=sum(ins) / len(ins), asp73_mean=stt.mean(dasp),
                                asp73_contact_fraction=sum(x <= 3.5 for x in dasp) / len(dasp),
                                centroid_mean=stt.mean(dcen), centroid_min=min(dcen), centroid_max=max(dcen),
                                loop_min=min(dloop) if dloop else None, loop_max=max(dloop) if dloop else None)

# ---------------------------------------------------------------------------------------- redock RMSD
for key, per in JOBS.items():
    if '/rli_anion__' not in key:
        continue
    r_ = [rmsd_sym(m[0][1], REF) for m in per]
    vals = [x[0] for x in r_]
    out['rmsd'][key] = dict(ref=f'dock/{C1}/ref_RLI.sdf', n=len(vals), mean=stt.mean(vals), sd=stt.stdev(vals),
                            min=min(vals), max=max(vals), f_le_1=sum(v <= 1.0 for v in vals) / len(vals),
                            f_le_2=sum(v <= 2.0 for v in vals) / len(vals), n_mappings=r_[0][1], per_seed=vals)


# ---------------------------------------------------------------------------------------- clorobiocin reference
def chain_a(pdb_id):
    with gzip.open(DOCK / 'pdb' / f'{pdb_id}.cif.gz', 'rt') as fh:
        st = gemmi.make_structure_from_block(gemmi.cif.read_string(fh.read()).sole_block())
    ca, het = [], {}
    for ch in st[0]:
        if ch.name != 'A':
            continue
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info and info.is_amino_acid():
                a = res.find_atom('CA', '*')
                if a:
                    ca.append((res.name, np.array([a.pos.x, a.pos.y, a.pos.z]), res.seqid.num))
            elif res.name != 'HOH':
                het[res.name] = [(a.element.name, np.array([a.pos.x, a.pos.y, a.pos.z])) for a in res
                                 if a.element.name != 'H']
    return ca, het


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


ca_k, het_k = chain_a('1KZN')
ca_d, _ = chain_a('4DUH')
sm = difflib.SequenceMatcher(None, [c[0] for c in ca_k], [c[0] for c in ca_d], autojunk=False)
pairs = [(i + k, j + k) for i, j, n in sm.get_matching_blocks() for k in range(n)]
P = np.array([ca_k[i][1] for i, _ in pairs])
Q = np.array([ca_d[j][1] for _, j in pairs])
R, t = kabsch(P, Q)
fit = (R @ P.T).T + t
ca_rmsd = float(np.sqrt(((fit - Q) ** 2).sum(1).mean()))
cbn_own = [(el, R @ x + t) for el, x in het_k['CBN']]
cbn_file = sdf_atoms(DOCK / C1 / 'ref_cross_X_1KZN_CBN.sdf')
dev = rmsd_sym(cbn_own, cbn_file)[0]
rec = receptor(C1)
RX = np.array([r_[3] for r_ in rec])
CX = np.array([a[1] for a in cbn_file])
D = np.linalg.norm(RX[:, None, :] - CX[None, :, :], axis=2)
i_min = int(D.min(1).argmin())
centre = np.array(json.loads((DOCK / C1 / 'box.json').read_text())['center'])
# residues of 4DUH (deposited chain A) that the alignment leaves without a partner in the 1KZN model
paired_d = {ca_d[j][2] for _, j in pairs}
unpaired_d = sorted(c[2] for c in ca_d if c[2] not in paired_d)
hit = [rec[a] for a in np.where(D < 2.5)[0]]                 # one entry per pair below 2.5 A
by_res = {}
for r_ in hit:
    by_res[f'{r_[1]}{r_[0]}'] = by_res.get(f'{r_[1]}{r_[0]}', 0) + 1
out['clorobiocin'] = dict(n_ca_pairs=len(pairs), ca_rmsd=ca_rmsd, own_vs_file_rmsd=dev,
                          min_heavy_heavy=float(D.min()),
                          closest=f'{rec[i_min][1]}{rec[i_min][0]}:{rec[i_min][2]}',
                          n_pairs_below_2p5=int((D < 2.5).sum()), n_pairs_below_3p0=int((D < 3.0).sum()),
                          pairs_below_2p5_by_residue=by_res,
                          residues_of_4DUH_without_counterpart_in_1KZN=unpaired_d,
                          n_pairs_below_2p5_without_counterpart=sum(1 for r_ in hit if r_[0] in unpaired_d),
                          residues_below_2p5_without_counterpart=sorted({f'{r_[1]}{r_[0]}' for r_ in hit if r_[0] in unpaired_d},
                                                                        key=lambda x: int(x[3:])),
                          margin_20=float(10.0 - np.abs(CX - centre).max()), n_heavy=len(cbn_file))

# ---------------------------------------------------------------------------------------- loop builds
out['loop'] = {}
for target in (C1, C1B):
    rec = receptor(target)
    L = np.array([r_[3] for r_ in rec if 109 <= r_[0] <= 116])
    out['loop'][target] = dict(n_heavy=len(L), to_inhibitor=float(np.linalg.norm(L[:, None] - REFX[None], axis=2).min()),
                               in_box=int((np.abs(L - centre) <= 10.0).all(axis=1).sum()))
a_, b_ = receptor(C1), {(r_[0], r_[2]): r_[3] for r_ in receptor(C1B)}
dl = [dist(r_[3], b_[(r_[0], r_[2])]) for r_ in a_ if 109 <= r_[0] <= 116]
do = [dist(r_[3], b_[(r_[0], r_[2])]) for r_ in a_ if not 109 <= r_[0] <= 116]
out['loop']['build_difference'] = dict(loop_rms=float(np.sqrt(np.mean(np.square(dl)))), loop_max=max(dl),
                                       other_rms=float(np.sqrt(np.mean(np.square(do)))), other_max=max(do))
out['loop']['closest_pose'] = min(v['loop_min'] for k, v in out['site'].items() if k.startswith(C1 + '/'))

# ---------------------------------------------------------------------------------------- size relation, targets
D9 = json.loads((OUT / 'd09_dock.json').read_text())['scores']
out['size'], out['targets'] = {}, {}
for tkey, rec_, sc in (('vina', C1, 'vina'), ('vinardo', C1, 'vinardo'), ('vina_water', C1W, 'vina')):
    rows = [(k.split('/runs/')[1].split('__')[0], v) for k, v in out['scores'].items()
            if k.startswith(rec_ + '/') and k.split('/runs/')[1].split('__')[1] == sc]
    if len(rows) > 2:
        ha = np.array([v['heavy_atoms'] for _, v in rows], float)
        sm_ = np.array([v['mean'] for _, v in rows])
        slope, icpt = np.polyfit(ha, sm_, 1)
        out['size'][tkey] = dict(n=len(rows), r_score_vs_HA=float(np.corrcoef(ha, sm_)[0, 1]),
                                 r_LE_vs_HA=float(np.corrcoef(ha, -sm_ / ha)[0, 1]),
                                 slope=float(slope), intercept=float(icpt),
                                 order=[k for k, _ in sorted(rows, key=lambda kv: kv[1]['mean'])])
    c1 = dict(rows)
    # like for like: the eight species of the production matrix, in the control receptor and in each structure
    tgs = (('T1_3TZF', '3TZF'), ('T5_7S2K', '7S2K'), ('T3_6RL9', '6RL9'), ('T4_3HS4', '3HS4'))
    com = [sp for sp in SPECIES8 if sp in c1 and all(f'{tg}/{sp}__{sc}__prod' in D9 for tg, _ in tgs)]
    if len(com) >= 4:
        ha8 = np.array([c1[sp]['heavy_atoms'] for sp in com], float)
        out['size'].setdefault(tkey, {})['common'] = dict(
            n=len(com), species=com,
            r=dict({'4DUH': float(np.corrcoef(ha8, [c1[sp]['mean'] for sp in com])[0, 1])},
                   **{pdb: float(np.corrcoef(ha8, [D9[f'{tg}/{sp}__{sc}__prod']['mean'] for sp in com])[0, 1])
                      for tg, pdb in tgs}))
    for tg, pdb in (('T1_3TZF', '3TZF'), ('T5_7S2K', '7S2K'), ('T3_6RL9', '6RL9'), ('T4_3HS4', '3HS4')):
        both = [(sp, c1[sp]['mean'], D9[f'{tg}/{sp}__{sc}__prod']['mean']) for sp in SPECIES8
                if sp in c1 and f'{tg}/{sp}__{sc}__prod' in D9]
        if len(both) >= 4:
            d_ = np.array([b[1] - b[2] for b in both])
            out['targets'][f'{tkey}/{pdb}'] = dict(n=len(both), mean_difference=float(d_.mean()),
                                                   min_difference=float(d_.min()), max_difference=float(d_.max()),
                                                   pearson_r=float(np.corrcoef([b[1] for b in both], [b[2] for b in both])[0, 1]),
                                                   per_species={b[0]: float(b[1] - b[2]) for b in both})

# neutral sulfanilamide, Vina: every variant of the control receptor against every production structure
sv = {t: out['scores'][f'{t}/runs/sfa_neutral__vina__prod']['mean'] for t in (C1, C1B, C1N, C1W)
      if f'{t}/runs/sfa_neutral__vina__prod' in out['scores']}
st_ = {pdb: D9[f'{tg}/sfa_neutral__vina__prod']['mean'] for tg, pdb in (('T1_3TZF', '3TZF'), ('T5_7S2K', '7S2K'),
                                                                         ('T3_6RL9', '6RL9'), ('T4_3HS4', '3HS4'))}
out['sfa_variants'] = dict(control=sv, production=st_, control_min=min(sv.values()), control_max=max(sv.values()),
                           max_abs_difference=max(abs(a - b) for a in sv.values() for b in st_.values()))

# ---------------------------------------------------------------------------------------- whole-fragment search
for d in sorted((DOCK / C1 / 'runs').glob('*__whole')):
    if not (d / 'results.json').exists():
        continue
    per = run(C1, d.name)
    top_in, any_in, gap = [], [], []
    for m in per:
        ins = [float(np.linalg.norm(REFX - np.array([a[1] for a in mm[1]]).mean(0), axis=1).min()) <= SITE_CUT for mm in m]
        top_in.append(ins[0])
        any_in.append(any(ins))
        if any(ins):
            gap.append(m[ins.index(True)][0] - m[0][0])
    w = dict(score_stats(per), n_top_in_site=int(sum(top_in)), n_any_in_site=int(sum(any_in)),
             max_gap_in_site_mode=max(gap) if gap else None)
    if d.name.startswith('rli_anion'):
        vals = [rmsd_sym(m[0][1], REF)[0] for m in per]
        w.update(rmsd_mean=stt.mean(vals), rmsd_sd=stt.stdev(vals) if len(vals) > 1 else 0.0, rmsd_max=max(vals),
                 f_le_1=sum(v <= 1.0 for v in vals) / len(vals))
    out['whole'][f'{C1}/runs/{d.name}'] = w

# ---------------------------------------------------------------------------------------- against the docking analysis
CS = json.loads((DOCK / 'control' / 'control_summary.json').read_text())
cmp_s = [abs(v['mean'] - CS['jobs'][k]['score']['mean']) for k, v in out['scores'].items() if k in CS['jobs']]
cmp_sd = [abs(v['sd'] - CS['jobs'][k]['score']['sd']) for k, v in out['scores'].items() if k in CS['jobs']]
cmp_r = [abs(v['mean'] - CS['validation'][k]['stats']['rmsd_top_mean']) for k, v in out['rmsd'].items()
         if k in CS['validation'] and CS['validation'][k]['stats'].get('rmsd_available')]
out['max_abs_score_mean_diff_vs_summary'] = max(cmp_s)
out['max_abs_score_sd_diff_vs_summary'] = max(cmp_sd)
out['max_abs_rmsd_mean_diff_vs_summary'] = max(cmp_r) if cmp_r else None
out['n_jobs_compared'] = len(cmp_s)
dump('d15_control.json', out)
print(f"jobs {len(cmp_s)}; max |d score mean| {max(cmp_s):.2e}; max |d SD| {max(cmp_sd):.2e}; "
      f"max |d RMSD| {(max(cmp_r) if cmp_r else float('nan')):.2e}; CBN transfer own-vs-file {dev:.3f} A "
      f"({len(pairs)} CA, {ca_rmsd:.3f} A)")
