"""Gate G-R, step 1: the geometries of the test set (md/PLAN.md, line of 6 Oct 2026, 3:02 AM IST, point 2).

    /root/cqc-sfa-tools/md/bin/python md/qm/gr_points.py            (NumPy only; no random numbers; seconds)

Writes md/qm/gr/points.json, points.csv, scan_ns_contacts.csv, scan_znn_contacts.csv and geom/<id>.xyz,
geom/<id>.pdb.  A second run gives the same files (the script compares and says so).  A geometry file that is
listed in md/qm/gr/SHA256SUMS must come out with the checksum listed there, or the script stops: quantum-
chemical energies exist for these files.

THE MODEL.  The 64 atoms of the small model of MCPB.py (md/qm/mcpb/CA2_small.pdb, frame of 6RL9): three
ethyl-imidazoles (45 atoms), the zinc, the sulfonamidate anion (18 atoms).  Zinc and imidazoles are FIXED in
every point; only the ligand changes.  Thr199 OG1 is not in the model; its position (hydrogen-bond acceptor of
N2-H) is read from the structure and used to place the hydrogen and to report the hydrogen bond.

NAMES.  z = unit vector N2->Zn, u = N2->S, h = N2->H.
  phi   out-of-plane angle of the N2-S bond: the angle between u and the plane through z and h, positive on
        the side of z x h.  +49.35 deg in start X, 0 for a planar nitrogen, negative in the configuration of R.
  psi   direction of the projection of u in that plane, counted from z towards h.
  tau   torsion Zn-N2-S-C4 (the turn about N2-S as seen from the zinc).
  B     the rigid body: the 17 ligand atoms other than HN21, with the internal geometry of the crystal pose
        (heavy atoms of SAN 304 of 6RL9, ring and aniline hydrogens as built for start X).
  T     the construction of the hydrogen used for all three starts (md/SYSTEM.md, section 5): N2 tetrahedral
        with S and Zn as two substituents, N-H 1.021 A, of the two free directions the one nearer to OG1.

THE POINTS
 (b) X      start X: the geometry MCPB.py wrote (identical to md/system/ligand_X.pdb).
     R, N   starts R and N as built (md/system/ligand_R.pdb, ligand_N.pdb): docked poses of the neutral
            molecule, so S-N 1.677 A and Zn-N2 2.21-2.24 A.
     Rx, Nx "R-like", "N-like": the body B fitted on the 11 heavy atoms of the start (least squares, rotation
            and translation), then moved along the line Zn->N2 until Zn-N2 has the crystal length; hydrogen
            by T.  Only the rigid-body placement of the crystal ligand differs from X.  The same recipe
            applied to start X returns X itself (checked), so the "X-like" point is X.
 (a) a1..a7 inversion at N2 from X (a1) to Rx (a7).  Zn-N2, N2-S, N2-H and B are kept.  One parameter s:
              z(s), h(s)  great-circle interpolation between their directions in X and in Rx
                          (N2 moves on the sphere about the zinc; H keeps pointing to OG1),
              phi(s), psi(s), tau(s)  linear between their values in X and in Rx,
            and the seven points are phi = phi_X * (1, 2/3, 1/3), 0, phi_Rx * (1/3, 2/3, 1), with
            s = (phi_X - phi) / (phi_X - phi_Rx).  a4 has a planar nitrogen (angle sum 360).  A path with N2
            fixed at its place in X cannot end at R: with the turn about N2-S of R the ring would then stand
            in the imidazole of His96 (C...C 1.9 A); R has N2 moved by 0.25 A on the sphere about the zinc.
 (c) c1_*   turn about N2-S in the configuration of X: B without N2 turned rigidly about the axis N2->S by
            delta, everything else as in X (tau = tau_X + delta).
     c2_*   the same from Rx.
            A turn is "possible" when, in the model, no ligand atom is closer to an imidazole atom than
            2.8 A (heavy-heavy), 2.1 A (heavy-hydrogen) or 1.9 A (hydrogen-hydrogen).  Taken are the five
            possible multiples of 30 deg nearest to the reference (ties: the negative delta first); where
            fewer than five are possible, possible multiples of 15 deg are added, the nearest to the
            reference first.  scan_ns_contacts.csv holds the whole turn in 5 deg steps.
 (d) d1..d6 hydrogen moved, heavy atoms as in X.  Not asked for by the plan; added because in (a) to (c) the
            angle Zn-N2-H hardly changes (105.8-110.4 deg), so that a fit could not tell its term from a
            constant.  omega = out-of-plane angle of N2-H against the plane Zn-N2-S (54.7 deg in X).
              d1      omega = 0: H in the plane Zn-N2-S (planar nitrogen by moving H; both H angles 121.2)
              d2      omega = -omega_X: H on the other side of that plane (all three angles as in X)
              d3, d4  H turned about the axis N2->S until Zn-N2-H is 12 deg larger, smaller (S-N2-H kept)
              d5, d6  H turned about the axis N2->Zn until S-N2-H is 12 deg larger, smaller (Zn-N2-H kept)
            d3 to d6 move ONE of the two hydrogen angles each: moves that change both together (as d1, or a
            rock in the bisecting plane) give only the sum of the two force constants.  In d1 and d2 the
            hydrogen no longer points to OG1; these geometries test the terms, not the pocket.

CONTACTS.  For every point: the closest ligand...imidazole distances in the model (heavy-heavy,
heavy-hydrogen, hydrogen-hydrogen), Zn...O1/O2, and the closest distance to a protein atom that is NOT in the
model (md/qm/mcpb/CA2_H.pdb, hydrogens as tleap placed them).  Residues carry the numbers of 6RL9.  Flags:
  model_contact    below the three limits above (such a grid point is not chosen in (c); elsewhere reported)
  protein_contact  a ligand heavy atom within 3.0 A of a protein heavy atom outside the model when a carbon
                   or the sulfur is involved, or within 2.4 A when both are N or O (hydrogen bonds are
                   2.5-3.0 A): the rigid protein of 6RL9 would not allow the geometry.  Contacts with protein
                   hydrogens are listed and not used (the hydrogens are tleap's, not relaxed).  Not dropped:
                   the model does not contain these atoms, so the comparison of the force field with quantum
                   chemistry is valid there; the tables give the gate numbers with and without them.
  drop             a heavy-heavy distance below 2.0 A in the model: not computed.

ADDED on 6 Oct 2026 after the review of the gate (md/PLAN.md, entry of 8:56 AM IST, point 2; review finding 6).
None of them changes a geometry of the first 26 points (the script checks the files against SHA256SUMS).
 (e) e_*    the whole ligand of X, hydrogen included, turned rigidly about the axis Zn->N2 by delta (sense: that
            of the torsion NE2(His94)-Zn-N2-S).  Bond lengths, the three angles at N2, the N(His)-Zn-N2 angles
            and every torsion inside the ligand stay; only the torsions through the zinc change, which the
            force field sets to zero, and the contacts.  The hydrogen bond to OG1 is given up.  Taken: the
            possible multiples of 15 deg with |delta| <= 30, and for +120 and -120 (the two other staggered
            places against the three histidines) that turn if possible, else the nearest possible multiple
            of 15 deg within 30 deg of it.  scan_znn_contacts.csv holds the whole turn in 5 deg steps.
 (f) f_p08, f_m08   the angle Zn-N2-S of X made 8 deg larger, smaller: B without N2 turned about the normal of
            the plane Zn-N2-S through N2; the hydrogen put where Zn-N2-H and S-N2-H have their values of X, on
            the same side of that plane.  Of the three angles at N2 only Zn-N2-S changes; Zn-N2-S-C4 is kept.
 (g) g1..g6 the six hydrogen moves of (d), made on Rx instead of X (the other configuration at N2).

ORDER: X, R, Rx, the planar point a4, then a3, a5, a2, a6, then d1 and d3 to d6 (without them a fit cannot
place the two hydrogen angles), d2, Nx, the (c) points with the smallest turns first and the two
configurations in turn, N; then the added points (e), (f), (g).  (Until 6 Oct 2026, 9 AM IST, d2 and Nx came
after the turns; review finding 11.)

SETS (md/PLAN.md, entry of 6 Oct 2026, 8:56 AM IST, point 2; used by md/qm/gr_fit.py).
  part_A  the terms at N2: X, the inner points of the path a2 to a6, (d), (e), (f), and Rx with (g).
  part_B  contacts: X, Rx, Nx and the turns (c).  Reported with the split by energy part, never fitted.
  base    the point a geometry was made from by moving atoms at N2 or by a turn: X, or Rx for (g) and c2_*.
  report_only   R and N as built (docked geometries of the neutral molecule).
  gate_plan, gate_rigid, fit: the sets of the first version of the gate (entries of 3:02 and 8:01 AM), kept
  so that the earlier numbers can be reproduced.
"""
import hashlib
import json
import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
MD = HERE.parent
OUT = HERE / 'gr'
SMALL = HERE / 'mcpb' / 'CA2_small.pdb'
PROT_H = HERE / 'mcpb' / 'CA2_H.pdb'
NH = 1.021
LIMITS = dict(heavy_heavy=2.8, heavy_h=2.1, h_h=1.9)
PROT_LIMITS = dict(polar_polar=2.4, other_heavy=3.0)
DROP_HEAVY = 2.0
GRID = 30


def read_pdb(path):
    out = []
    for l in pathlib.Path(path).read_text().splitlines():
        if l.startswith(('ATOM', 'HETATM')):
            name = l[12:16].strip()
            resn = l[17:20].strip()
            el = 'Zn' if resn == 'ZN' else name[0]
            out.append(dict(serial=int(l[6:11]), name=name, resn=resn, resi=int(l[22:26]), el=el,
                            xyz=np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return out


def unit(v):
    return v / np.linalg.norm(v)


def angle(a, b, c):
    return float(np.degrees(np.arccos(np.clip(np.dot(unit(a - b), unit(c - b)), -1, 1))))


def dihedral(a, b, c, d):
    b1 = unit(c - b)
    v = (a - b) - np.dot(a - b, b1) * b1
    w = (d - c) - np.dot(d - c, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def rot(axis, deg):
    a = unit(axis)
    t = np.radians(deg)
    k = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(t) * k + (1 - np.cos(t)) * (k @ k)


def slerp(a, b, s):
    w = np.arccos(np.clip(np.dot(a, b), -1, 1))
    if w < 1e-12:
        return a.copy()
    return (np.sin((1 - s) * w) * a + np.sin(s * w) * b) / np.sin(w)


def kabsch(p, q):
    """Rotation r and the two centroids: r (p - pc) + qc is the least-squares fit of p on q."""
    pc, qc = p.mean(0), q.mean(0)
    u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    return vt.T @ np.diag([1.0, 1.0, d]) @ u.T, pc, qc


def wrap(x):
    return (x + 180.0) % 360.0 - 180.0


def deposited_numbers():
    """tleap residue index -> residue number of 6RL9 (md/system/residue_map.csv)."""
    out = {}
    rows = (MD / 'system' / 'residue_map.csv').read_text().splitlines()
    head = rows[0].split(',')
    for r in rows[1:]:
        f = dict(zip(head, r.split(',')))
        out[int(f['tleap_index'])] = f"{f['deposited_name'].split()[0].capitalize()}{f['deposited_number']}"
    return out


class Model:
    def __init__(self):
        sm = read_pdb(SMALL)
        dep_no = deposited_numbers()
        assert len(sm) == 64
        self.atoms = sm
        self.lig = [k for k, a in enumerate(sm) if a['resn'] == 'SAN']
        self.fix = [k for k, a in enumerate(sm) if a['resn'] != 'SAN']
        assert len(self.lig) == 18 and self.lig == list(range(46, 64))
        self.names = [sm[k]['name'] for k in self.lig]
        self.el = [sm[k]['el'] for k in self.lig]
        self.x0 = np.array([sm[k]['xyz'] for k in self.lig])
        self.fxyz = np.array([sm[k]['xyz'] for k in self.fix])
        self.fel = [sm[k]['el'] for k in self.fix]
        self.flab = [f"{dep_no.get(sm[k]['resi'], sm[k]['resn'])}:{sm[k]['name']}" for k in self.fix]
        self.zn = np.array([a['xyz'] for a in sm if a['el'] == 'Zn'][0])
        # the methyl groups that stand for the backbone (atoms H1, CH3, H2, H3 of each ethyl-imidazole)
        self.caps = [n for n, k in enumerate(self.fix) if sm[k]['name'] in ('H1', 'CH3', 'H2', 'H3')]
        assert len(self.caps) == 12
        ph = read_pdb(PROT_H)
        og = [a for a in ph if a['resn'] == 'THR' and a['resi'] == 197 and a['name'] == 'OG1']
        assert len(og) == 1
        self.og1 = og[0]['xyz']
        # the same atom in the deposited numbering, as a check of the residue number
        dep = [a for a in read_pdb(MD / 'system' / 'protein_md.pdb')
               if a['resn'] == 'THR' and a['resi'] == 199 and a['name'] == 'OG1']
        assert len(dep) == 1 and np.abs(dep[0]['xyz'] - self.og1).max() < 1e-9
        in_model = {a['serial'] for a in sm}
        prot = [a for a in ph if a['serial'] not in in_model and a['resn'] not in ('ZN', 'SAN')]
        self.pxyz = np.array([a['xyz'] for a in prot])
        self.pel = [a['el'] for a in prot]
        self.plab = [f"{dep_no[a['resi']]}:{a['name']}" for a in prot]
        self.i = {n: k for k, n in enumerate(self.names)}
        self.body = [k for k, n in enumerate(self.names) if n != 'HN21']
        self.turn = [k for k, n in enumerate(self.names) if n not in ('HN21', 'N2')]
        self.heavy = [k for k, e in enumerate(self.el) if e != 'H']
        self.d_zn_n = float(np.linalg.norm(self.x0[self.i['N2']] - self.zn))
        self.n94 = [a['xyz'] for a in sm if a['resi'] == 93 and a['name'] == 'NE2'][0]
        self.d_n_s = float(np.linalg.norm(self.x0[self.i['S']] - self.x0[self.i['N2']]))

    # ------------------------------------------------------------------------------------------ constructions
    def build_h(self, n2, s):
        a, b = unit(s - n2), unit(self.zn - n2)
        bis = unit(-(a + b))
        perp = unit(np.cross(a, b))
        half = np.radians(109.47 / 2.0)
        cand = [n2 + NH * (np.cos(half) * bis + sg * np.sin(half) * perp) for sg in (1, -1)]
        return cand[int(np.argmin([np.linalg.norm(c - self.og1) for c in cand]))]

    def start(self, tag):
        lig = {a['name']: a['xyz'] for a in read_pdb(MD / 'system' / f'ligand_{tag}.pdb')}
        return np.array([lig[n] for n in self.names])

    def like(self, built):
        r, pc, qc = kabsch(self.x0[self.heavy], built[self.heavy])
        y = (r @ (self.x0 - pc).T).T + qc
        fit_rmsd = float(np.sqrt(((y[self.heavy] - built[self.heavy]) ** 2).sum(1).mean()))
        v = y[self.i['N2']] - self.zn
        y = y + (self.d_zn_n - np.linalg.norm(v)) * unit(v)
        y[self.i['HN21']] = self.build_h(y[self.i['N2']], y[self.i['S']])
        return y, fit_rmsd

    def frame(self, c):
        n2 = c[self.i['N2']]
        z, h, u = unit(self.zn - n2), unit(c[self.i['HN21']] - n2), unit(c[self.i['S']] - n2)
        eb = unit(h - np.dot(h, z) * z)
        n = np.cross(z, eb)
        phi = np.degrees(np.arcsin(np.clip(np.dot(u, n), -1, 1)))
        p = unit(u - np.dot(u, n) * n)
        psi = np.degrees(np.arctan2(np.dot(p, eb), np.dot(p, z)))
        tau = dihedral(self.zn, n2, c[self.i['S']], c[self.i['C4']])
        return dict(z=z, h=h, phi=float(phi), psi=float(psi), tau=tau)

    def place(self, z, h, phi, psi, tau):
        """The ligand with N2 at Zn - d z, H at N2 + 1.021 h, S in the direction (phi, psi), turn tau."""
        n2 = self.zn - self.d_zn_n * z
        eb = unit(h - np.dot(h, z) * z)
        n = np.cross(z, eb)
        p = np.cos(np.radians(psi)) * z + np.sin(np.radians(psi)) * eb
        u = np.cos(np.radians(phi)) * p + np.sin(np.radians(phi)) * n
        # body: N2 on its place, the N2->S direction of the crystal body turned on u, then the turn about u
        b0 = self.x0 - self.x0[self.i['N2']]
        u0 = unit(b0[self.i['S']])
        ax = np.cross(u0, u)
        if np.linalg.norm(ax) < 1e-12:
            r1 = np.eye(3)
        else:
            r1 = rot(ax, np.degrees(np.arctan2(np.linalg.norm(ax), np.dot(u0, u))))
        c = (r1 @ b0.T).T + n2
        cur = dihedral(self.zn, n2, c[self.i['S']], c[self.i['C4']])
        c = self.turn_ns(c, wrap(tau - cur))
        c[self.i['HN21']] = n2 + NH * h
        assert abs(wrap(dihedral(self.zn, n2, c[self.i['S']], c[self.i['C4']]) - tau)) < 1e-7
        return c

    def turn_ns(self, c, delta):
        """B without N2 turned about N2->S; the sense is that of the torsion Zn-N2-S-C4."""
        n2, s = c[self.i['N2']], c[self.i['S']]
        t0 = dihedral(self.zn, n2, s, c[self.i['C4']])
        for sg in (1, -1):
            d = c.copy()
            d[self.turn] = (rot(s - n2, sg * delta) @ (c[self.turn] - n2).T).T + n2
            if abs(wrap(dihedral(self.zn, n2, d[self.i['S']], d[self.i['C4']]) - t0 - delta)) < 1e-7:
                return d
        raise RuntimeError('turn about N2-S failed')

    def move_h(self, c, omega, rho):
        n2 = c[self.i['N2']]
        z, u = unit(self.zn - n2), unit(c[self.i['S']] - n2)
        nzs = unit(np.cross(u, z))
        q = rot(nzs, rho) @ unit(-(z + u))
        d = c.copy()
        d[self.i['HN21']] = n2 + NH * (np.cos(np.radians(omega)) * q + np.sin(np.radians(omega)) * nzs)
        return d

    def turn_h(self, c, which, delta):
        """H turned about N2->S (which = 'Zn_N2_H': S-N2-H is kept) or about N2->Zn (which = 'S_N2_H':
        Zn-N2-H is kept) until the named angle has changed by delta degrees; the smaller turn of the two."""
        n2, hn = c[self.i['N2']], c[self.i['HN21']]
        other = self.zn if which == 'Zn_N2_H' else c[self.i['S']]
        axis = (c[self.i['S']] if which == 'Zn_N2_H' else self.zn) - n2
        base = angle(other, n2, hn)

        def f(alpha):
            return angle(other, n2, n2 + rot(axis, alpha) @ (hn - n2)) - base - delta
        best = None
        for sg in (1.0, -1.0):
            lo, hi = 0.0, 0.0
            while abs(hi) < 90.0 and f(lo) * f(hi) > 0:
                hi += sg * 1.0
            if f(lo) * f(hi) > 0:
                continue
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if f(lo) * f(mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            alpha = 0.5 * (lo + hi)
            if best is None or abs(alpha) < abs(best):
                best = alpha
        if best is None:
            raise RuntimeError(f'no turn of H gives {which} {delta:+.1f} deg')
        d = c.copy()
        d[self.i['HN21']] = n2 + rot(axis, best) @ (hn - n2)
        return d, best

    def turn_znn(self, c, delta):
        """The whole ligand turned rigidly about Zn->N2; the sense is that of the torsion NE2(His94)-Zn-N2-S."""
        n2 = c[self.i['N2']]
        t0 = dihedral(self.n94, self.zn, n2, c[self.i['S']])
        for sg in (1, -1):
            d = (rot(n2 - self.zn, sg * delta) @ (c - n2).T).T + n2
            if abs(wrap(dihedral(self.n94, self.zn, n2, d[self.i['S']]) - t0 - delta)) < 1e-7:
                return d
        raise RuntimeError('turn about Zn-N2 failed')

    def bend_zns(self, c, delta):
        """Zn-N2-S changed by delta in its plane (B without N2 turned about the normal of the plane through
        N2); the hydrogen where Zn-N2-H and S-N2-H keep their values, on the same side of the plane."""
        n2, s, hn = c[self.i['N2']], c[self.i['S']], c[self.i['HN21']]
        z, u, h = unit(self.zn - n2), unit(s - n2), unit(hn - n2)
        nrm = unit(np.cross(z, u))
        a0 = angle(self.zn, n2, s)
        for sg in (1, -1):
            d = c.copy()
            d[self.turn] = (rot(nrm, sg * delta) @ (c[self.turn] - n2).T).T + n2
            if abs(angle(self.zn, n2, d[self.i['S']]) - a0 - delta) < 1e-7:
                break
        else:
            raise RuntimeError('bend of Zn-N2-S failed')
        u2 = unit(d[self.i['S']] - n2)
        a1, a2, g = float(np.dot(h, z)), float(np.dot(h, u)), float(np.dot(z, u2))
        p, q = (a1 - a2 * g) / (1 - g * g), (a2 - a1 * g) / (1 - g * g)
        flat = p * z + q * u2
        out2 = 1.0 - float(np.dot(flat, flat))
        if out2 <= 0:
            raise RuntimeError('no place for the hydrogen with both angles kept')
        n2u = unit(np.cross(z, u2))
        assert np.dot(n2u, nrm) > 0.999999
        d[self.i['HN21']] = n2 + np.linalg.norm(hn - n2) * (flat + np.sign(np.dot(h, nrm)) * np.sqrt(out2) * n2u)
        return d

    def omega_rho(self, c):
        n2 = c[self.i['N2']]
        z, u, h = unit(self.zn - n2), unit(c[self.i['S']] - n2), unit(c[self.i['HN21']] - n2)
        nzs = unit(np.cross(u, z))
        om = np.degrees(np.arcsin(np.clip(np.dot(h, nzs), -1, 1)))
        q = unit(h - np.dot(h, nzs) * nzs)
        q0 = unit(-(z + u))
        rho = np.degrees(np.arctan2(np.dot(np.cross(q0, q), nzs), np.dot(q0, q)))
        return float(om), float(rho)

    # --------------------------------------------------------------------------------------------- description
    def describe(self, c):
        i, zn = self.i, self.zn
        n2, s, hn, c4 = c[i['N2']], c[i['S']], c[i['HN21']], c[i['C4']]
        f = self.frame(c)
        om, rho = self.omega_rho(c)
        z, h, u = unit(zn - n2), unit(hn - n2), unit(s - n2)
        d = dict(
            Zn_N2=np.linalg.norm(zn - n2), N2_S=np.linalg.norm(s - n2), N2_H=np.linalg.norm(hn - n2),
            Zn_N2_S=angle(zn, n2, s), Zn_N2_H=angle(zn, n2, hn), S_N2_H=angle(s, n2, hn),
            signed_volume=float(np.dot(z, np.cross(h, u))), phi=f['phi'], psi=f['psi'], omega=om, rho=rho,
            Zn_N2_S_C4=f['tau'], H_N2_S_C4=dihedral(hn, n2, s, c4),
            H_N2_S_O1=dihedral(hn, n2, s, c[i['O1']]), H_N2_S_O2=dihedral(hn, n2, s, c[i['O2']]),
            Zn_N2_S_O1=dihedral(zn, n2, s, c[i['O1']]), Zn_N2_S_O2=dihedral(zn, n2, s, c[i['O2']]),
            C4_S_N2=angle(c4, s, n2), C3_C4_S_N2=dihedral(c[i['C3']], c4, s, n2),
            Zn_O1=np.linalg.norm(c[i['O1']] - zn), Zn_O2=np.linalg.norm(c[i['O2']] - zn),
            Zn_S=np.linalg.norm(s - zn),
            H_OG1=np.linalg.norm(hn - self.og1), N2_H_OG1=angle(n2, hn, self.og1),
            N2_OG1=np.linalg.norm(n2 - self.og1))
        d['angle_sum_N2'] = d['Zn_N2_S'] + d['Zn_N2_H'] + d['S_N2_H']
        d['N94_Zn_N2_S'] = dihedral(self.n94, zn, n2, s)
        for res, at, key in ((93, 'NE2', 'N94_Zn_N2'), (95, 'NE2', 'N96_Zn_N2'), (118, 'ND1', 'N119_Zn_N2')):
            p = [a['xyz'] for a in self.atoms if a['resi'] == res and a['name'] == at][0]
            d[key] = angle(p, zn, n2)
        d['rmsd_to_X_heavy'] = float(np.sqrt(((c[self.heavy] - self.x0[self.heavy]) ** 2).sum(1).mean()))
        d['closest_cap_atom'] = float(np.linalg.norm(c[:, None, :] - self.fxyz[self.caps][None, :, :], axis=2).min())
        return {k: round(float(v), 4) for k, v in d.items()}

    def contacts(self, c):
        best = dict(heavy_heavy=[9e9, ''], heavy_h=[9e9, ''], h_h=[9e9, ''])
        dm = np.linalg.norm(c[:, None, :] - self.fxyz[None, :, :], axis=2)
        for a in range(len(c)):
            for b in range(len(self.fxyz)):
                if self.fel[b] == 'Zn':
                    continue
                nh = (self.el[a] == 'H') + (self.fel[b] == 'H')
                t = best[('heavy_heavy', 'heavy_h', 'h_h')[nh]]
                if dm[a, b] < t[0]:
                    t[0], t[1] = float(dm[a, b]), f'{self.names[a]}...{self.flab[b]}'
        dp = np.linalg.norm(c[:, None, :] - self.pxyz[None, :, :], axis=2)
        pb = dict(polar_polar=[9e9, ''], other_heavy=[9e9, ''], any_h=[9e9, ''])
        for a in range(len(c)):
            for b in np.where(dp[a] < 6.0)[0]:
                if self.el[a] == 'H' or self.pel[b] == 'H':
                    key = 'any_h'
                elif self.el[a] in 'NO' and self.pel[b] in 'NO':
                    key = 'polar_polar'
                else:
                    key = 'other_heavy'
                if dp[a, b] < pb[key][0]:
                    pb[key][0], pb[key][1] = float(dp[a, b]), f'{self.names[a]}...{self.plab[b]}'
        out = dict(model={k: dict(d=round(v[0], 3), pair=v[1]) for k, v in best.items()},
                   protein_outside_model={k: (dict(d=round(v[0], 3), pair=v[1]) if v[0] < 9e8 else None)
                                          for k, v in pb.items()})
        flags = []
        if any(best[k][0] < LIMITS[k] for k in LIMITS):
            flags.append('model_contact')
        if any(pb[k][0] < PROT_LIMITS[k] for k in PROT_LIMITS):
            flags.append('protein_contact')
        if best['heavy_heavy'][0] < DROP_HEAVY:
            flags.append('drop')
        return out, flags


def choose_grid(scan):
    """Five possible turns: multiples of 30 deg nearest to the reference, filled up with multiples of 15 deg
    (delta = 0 excluded; ties: negative first)."""
    def possible(step):
        ok = [r['delta'] for r in scan
              if r['delta'] % step == 0 and r['delta'] != 0 and 'model_contact' not in r['flags']]
        return sorted(ok, key=lambda d: (abs(d), d))
    take = possible(GRID)[:5]
    for d in possible(GRID // 2):
        if len(take) >= 5:
            break
        if d not in take:
            take.append(d)
    return sorted(take)


def main():
    m = Model()
    OUT.mkdir(exist_ok=True)
    (OUT / 'geom').mkdir(exist_ok=True)
    pts = {}

    def add(pid, group, c, definition, **par):
        desc = m.describe(c)
        con, flags = m.contacts(c)
        pts[pid] = dict(id=pid, group=group, definition=definition, parameters=par, internal=desc, contacts=con,
                        flags=flags, coords=c)

    x = m.x0.copy()
    assert np.abs(x - m.start('X')).max() == 0.0, 'small model and md/system/ligand_X.pdb differ'
    add('X', 'b', x, 'start X, the geometry of the small model as MCPB.py wrote it')
    r_built, n_built = m.start('R'), m.start('N')
    add('R', 'b', r_built, 'start R as built (docked pose, bond lengths of the neutral molecule)')
    add('N', 'b', n_built, 'start N as built (docked pose, bond lengths of the neutral molecule)')
    xx, f0 = m.like(x)
    assert np.abs(xx - x).max() < 2e-3, 'the recipe of the like-points does not return X from X'
    rx, fr = m.like(r_built)
    nx, fn = m.like(n_built)
    add('Rx', 'b', rx, 'R-like: crystal ligand fitted on start R, Zn-N2 at the crystal length, H by T',
        fit_rmsd_heavy=round(fr, 4))
    add('Nx', 'b', nx, 'N-like: crystal ligand fitted on start N, Zn-N2 at the crystal length, H by T',
        fit_rmsd_heavy=round(fn, 4))
    like_check = dict(recipe_on_X_largest_shift_A=float(np.abs(xx - x).max()), fit_rmsd_on_X=f0)

    # ---- (a) the inversion path X -> Rx
    fx, fr_ = m.frame(x), m.frame(rx)
    # the parametrisation gives both ends back: the body to rounding, the hydrogen to the three decimals of
    # the files (N2-H is 1.021 A in the path; the coordinates of X are rounded to 0.001 A)
    for f, ref in ((fx, x), (fr_, rx)):
        back = m.place(f['z'], f['h'], f['phi'], f['psi'], f['tau'])
        assert np.abs(back[m.body] - ref[m.body]).max() < 1e-9 and np.abs(back - ref).max() < 2e-3
    dpsi, dtau = wrap(fr_['psi'] - fx['psi']), wrap(fr_['tau'] - fx['tau'])
    phis = [fx['phi'], fx['phi'] * 2 / 3, fx['phi'] / 3, 0.0, fr_['phi'] / 3, fr_['phi'] * 2 / 3, fr_['phi']]
    path = []
    for k, phi in enumerate(phis):
        s = (fx['phi'] - phi) / (fx['phi'] - fr_['phi'])
        c = m.place(slerp(fx['z'], fr_['z'], s), slerp(fx['h'], fr_['h'], s), phi, fx['psi'] + s * dpsi,
                    fx['tau'] + s * dtau)
        path.append(c)
        if k in (0, 6):
            ref = x if k == 0 else rx
            assert np.abs(c[m.body] - ref[m.body]).max() < 1e-9 and np.abs(c - ref).max() < 2e-3
            pts['X' if k == 0 else 'Rx']['parameters'].update(path_point=f'a{k + 1}', s=round(s, 6))
        else:
            add(f'a{k + 1}', 'a', c, 'inversion path X -> Rx' + (', planar nitrogen' if k == 3 else ''),
                s=round(s, 6), phi=round(phi, 4), psi=round(fx['psi'] + s * dpsi, 4),
                tau=round(wrap(fx['tau'] + s * dtau), 4))
    path_info = dict(
        phi_X=fx['phi'], phi_Rx=fr_['phi'], psi_X=fx['psi'], psi_Rx=fr_['psi'], tau_X=fx['tau'], tau_Rx=fr_['tau'],
        angle_between_z_X_and_z_Rx=float(np.degrees(np.arccos(np.dot(fx['z'], fr_['z'])))),
        angle_between_h_X_and_h_Rx=float(np.degrees(np.arccos(np.dot(fx['h'], fr_['h'])))),
        N2_shift_X_to_Rx_A=float(np.linalg.norm(x[m.i['N2']] - rx[m.i['N2']])))

    # why N2 has to move: N2 and H as in X, S on the other side (phi -> -phi), the turn about N2-S of start R
    tau_r = m.describe(r_built)['Zn_N2_S_C4']
    c_fixed = m.place(fx['z'], fx['h'], -fx['phi'], fx['psi'], tau_r)
    con_fixed, flags_fixed = m.contacts(c_fixed)
    n2_fixed = dict(what='N2 and H as in X, N2-S mirrored through the plane Zn-N2-H, Zn-N2-S-C4 of start R as built',
                    Zn_N2_S_C4=round(tau_r, 3), model_contacts=con_fixed['model'], flags=flags_fixed)

    # ---- (c) the turn about N2-S in both configurations
    scan_rows, chosen = [], {}
    ranges = {}
    for cfg, ref in (('c1', x), ('c2', rx)):
        scan = []
        for delta in range(-180, 180, 5):
            c = m.turn_ns(ref, delta)
            con, flags = m.contacts(c)
            d = m.describe(c)
            scan.append(dict(config=cfg, delta=delta, Zn_N2_S_C4=d['Zn_N2_S_C4'], H_N2_S_C4=d['H_N2_S_C4'],
                             heavy_heavy=con['model']['heavy_heavy']['d'], heavy_h=con['model']['heavy_h']['d'],
                             h_h=con['model']['h_h']['d'], pair_heavy_heavy=con['model']['heavy_heavy']['pair'],
                             Zn_O1=d['Zn_O1'], Zn_O2=d['Zn_O2'],
                             protein_polar_polar=(con['protein_outside_model']['polar_polar'] or {}).get('d'),
                             protein_other_heavy=(con['protein_outside_model']['other_heavy'] or {}).get('d'),
                             flags=flags))
        scan_rows += scan
        chosen[cfg] = choose_grid(scan)

        def runs(ok):
            out, cur = [], None
            for r in scan:
                if ok(r):
                    cur = [r['delta'], r['delta']] if cur is None else [cur[0], r['delta']]
                elif cur is not None:
                    out.append(cur)
                    cur = None
            if cur is not None:
                out.append(cur)
            return out
        ranges[cfg] = dict(
            possible_in_the_model=runs(lambda r: 'model_contact' not in r['flags']),
            also_without_protein_contact=runs(lambda r: not r['flags']),
            closest_heavy_heavy_in_the_model=min((dict(delta=r['delta'], d=r['heavy_heavy'], pair=r['pair_heavy_heavy'])
                                                  for r in scan), key=lambda v: v['d']))
        if cfg == 'c1':
            # the turn that would give the configuration of X the torsion Zn-N2-S-C4 of start R
            want = int(round(wrap(tau_r - fx['tau']) / 5.0)) * 5
            hit = [r for r in scan if r['delta'] == want][0]
            ranges[cfg]['turn_to_the_torsion_of_R'] = dict(delta=want, heavy_heavy=hit['heavy_heavy'],
                                                           pair=hit['pair_heavy_heavy'], flags=hit['flags'])
        for delta in chosen[cfg]:
            sign = 'm' if delta < 0 else 'p'
            add(f'{cfg}_{sign}{abs(delta):03d}', 'c', m.turn_ns(ref, delta),
                f"turn about N2-S by {delta:+d} deg from {'X' if cfg == 'c1' else 'Rx'}", delta=delta)

    # ---- (d) the hydrogen moved, heavy atoms as in X
    om_x, rho_x = m.omega_rho(x)
    assert om_x > 0 and abs(rho_x) < 0.05 and np.abs(m.move_h(x, om_x, 0.0) - x).max() < 2e-3
    for pid, om, text in (('d1', 0.0, 'H in the plane Zn-N2-S (planar nitrogen by moving H)'),
                          ('d2', -om_x, 'H on the other side of the plane Zn-N2-S')):
        add(pid, 'd', m.move_h(x, om, 0.0), text + '; heavy atoms as in X', omega=round(om, 4))
    for pid, which, delta in (('d3', 'Zn_N2_H', 12.0), ('d4', 'Zn_N2_H', -12.0), ('d5', 'S_N2_H', 12.0),
                              ('d6', 'S_N2_H', -12.0)):
        c, alpha = m.turn_h(x, which, delta)
        kept = 'S_N2_H' if which == 'Zn_N2_H' else 'Zn_N2_H'
        dsc = m.describe(c)
        assert abs(dsc[which] - pts['X']['internal'][which] - delta) < 1e-3
        assert abs(dsc[kept] - pts['X']['internal'][kept]) < 1e-3
        add(pid, 'd', c, f"H turned about N2->{'S' if which == 'Zn_N2_H' else 'Zn'} by {alpha:+.2f} deg: "
            f"{which.replace('_', '-')} {delta:+.0f} deg, {kept.replace('_', '-')} kept; heavy atoms as in X",
            angle=which, delta=delta, turn_deg=round(alpha, 4))

    # ---- added 6 Oct 2026: (e) turns about Zn-N2, (f) the bend of Zn-N2-S, (g) the hydrogen moves on Rx
    zscan = []
    for delta in range(-180, 180, 5):
        c = m.turn_znn(x, delta)
        con, flags = m.contacts(c)
        d = m.describe(c)
        zscan.append(dict(delta=delta, N94_Zn_N2_S=d['N94_Zn_N2_S'], heavy_heavy=con['model']['heavy_heavy']['d'],
                          pair_heavy_heavy=con['model']['heavy_heavy']['pair'], heavy_h=con['model']['heavy_h']['d'],
                          pair_heavy_h=con['model']['heavy_h']['pair'], h_h=con['model']['h_h']['d'],
                          pair_h_h=con['model']['h_h']['pair'], Zn_O1=d['Zn_O1'], Zn_O2=d['Zn_O2'], H_OG1=d['H_OG1'],
                          protein_polar_polar=(con['protein_outside_model']['polar_polar'] or {}).get('d'),
                          protein_other_heavy=(con['protein_outside_model']['other_heavy'] or {}).get('d'),
                          flags=flags))
    zok = sorted(r['delta'] for r in zscan if r['delta'] % 15 == 0 and r['delta'] != 0 and 'model_contact' not in r['flags'])
    zturns = [d for d in zok if abs(d) <= 30]
    for target in (120, -120):
        near = sorted((d for d in zok if abs(d - target) <= 30), key=lambda d: (abs(d - target), abs(d)))
        if near and near[0] not in zturns:
            zturns.append(near[0])
    zturns = sorted(zturns)
    e_ids = []
    for delta in zturns:
        pid = f"e_{'m' if delta < 0 else 'p'}{abs(delta):03d}"
        e_ids.append(pid)
        add(pid, 'e', m.turn_znn(x, delta), f'the ligand of X turned rigidly about Zn-N2 by {delta:+d} deg', delta=delta)
    e_ids = sorted(e_ids, key=lambda p: (abs(pts[p]['parameters']['delta']), pts[p]['parameters']['delta']))
    f_ids = []
    for delta in (8.0, -8.0):
        pid = f"f_{'m' if delta < 0 else 'p'}{abs(int(delta)):02d}"
        f_ids.append(pid)
        c = m.bend_zns(x, delta)
        dsc = m.describe(c)
        xi = pts['X']['internal']
        assert abs(dsc['Zn_N2_S'] - xi['Zn_N2_S'] - delta) < 1e-3 and abs(dsc['Zn_N2_H'] - xi['Zn_N2_H']) < 1e-3
        assert abs(dsc['S_N2_H'] - xi['S_N2_H']) < 1e-3 and abs(wrap(dsc['Zn_N2_S_C4'] - xi['Zn_N2_S_C4'])) < 1e-3
        assert dsc['signed_volume'] * xi['signed_volume'] > 0
        add(pid, 'f', c, f'Zn-N2-S {delta:+.0f} deg in its plane, Zn-N2-H and S-N2-H kept; from X', delta=delta)
    om_r, rho_r = m.omega_rho(rx)
    assert abs(rho_r) < 0.05 and np.abs(m.move_h(rx, om_r, 0.0) - rx).max() < 2e-3
    g_ids = ['g1', 'g2', 'g3', 'g4', 'g5', 'g6']
    for pid, om, text in (('g1', 0.0, 'H in the plane Zn-N2-S (planar nitrogen by moving H)'),
                          ('g2', -om_r, 'H on the other side of the plane Zn-N2-S')):
        add(pid, 'g', m.move_h(rx, om, 0.0), text + '; heavy atoms as in Rx', omega=round(om, 4))
    for pid, which, delta in (('g3', 'Zn_N2_H', 12.0), ('g4', 'Zn_N2_H', -12.0), ('g5', 'S_N2_H', 12.0),
                              ('g6', 'S_N2_H', -12.0)):
        c, alpha = m.turn_h(rx, which, delta)
        kept = 'S_N2_H' if which == 'Zn_N2_H' else 'Zn_N2_H'
        dsc = m.describe(c)
        assert abs(dsc[which] - pts['Rx']['internal'][which] - delta) < 1e-3
        assert abs(dsc[kept] - pts['Rx']['internal'][kept]) < 1e-3
        add(pid, 'g', c, f"H turned about N2->{'S' if which == 'Zn_N2_H' else 'Zn'} by {alpha:+.2f} deg: "
            f"{which.replace('_', '-')} {delta:+.0f} deg, {kept.replace('_', '-')} kept; heavy atoms as in Rx",
            angle=which, delta=delta, turn_deg=round(alpha, 4))

    c1, c2 = [f"c1_{'m' if d < 0 else 'p'}{abs(d):03d}" for d in chosen['c1']], \
             [f"c2_{'m' if d < 0 else 'p'}{abs(d):03d}" for d in chosen['c2']]

    def nearest_first(ids):
        return sorted(ids, key=lambda p: (abs(pts[p]['parameters']['delta']), pts[p]['parameters']['delta']))
    c1, c2 = nearest_first(c1), nearest_first(c2)
    mixed = [p for pair in zip(c1, c2) for p in pair] + c1[len(c2):] + c2[len(c1):]
    order = (['X', 'R', 'Rx', 'a4', 'a3', 'a5', 'a2', 'a6', 'd1', 'd3', 'd4', 'd5', 'd6', 'd2', 'Nx'] + mixed + ['N']
             + e_ids + f_ids + g_ids)
    assert sorted(order) == sorted(pts), (sorted(order), sorted(pts))
    order = [p for p in order if 'drop' not in pts[p]['flags']]

    sym = [a['el'] for a in m.atoms]
    resname = {('HID', 93): 'HD1', ('HID', 95): 'HD2', ('HIE', 118): 'HE1', ('ZN', 260): 'ZN1', ('SAN', 261): 'SN1'}
    listing = []
    # geometry files for which energies may exist: their text must not change (checksums of the first set)
    locked = {}
    if (OUT / 'SHA256SUMS').exists():
        for l in (OUT / 'SHA256SUMS').read_text().splitlines():
            if l.strip() and l.split()[1].startswith('gr/geom/') and l.split()[1].endswith('.xyz'):
                locked[l.split()[1][8:-4]] = l.split()[0]
    for n, pid in enumerate(order, 1):
        p = pts[pid]
        full = np.vstack([m.fxyz, p['coords']])
        assert len(full) == 64 and m.fix == list(range(46))
        xyz = [f'64', f"G-R point {pid} ({p['group']}): {p['definition']}; frame of 6RL9, Angstrom"]
        xyz += [f'{s:2s} {v[0]:14.8f} {v[1]:14.8f} {v[2]:14.8f}' for s, v in zip(sym, full)]
        text_xyz = '\n'.join(xyz) + '\n'
        if pid in locked and hashlib.sha256(text_xyz.encode()).hexdigest() != locked[pid]:
            raise SystemExit(f'geom/{pid}.xyz would change (checksum of md/qm/gr/SHA256SUMS): nothing written for it')
        (OUT / 'geom' / f'{pid}.xyz').write_text(text_xyz)
        pdb = [f"REMARK G-R point {pid}: {p['definition']}",
               'REMARK frame of 6RL9; atom order and names of the MM topology of the small model (md/qm/gr_mm.py),',
               'REMARK which is the order of the quantum-chemical input']
        for k, (a, v) in enumerate(zip(m.atoms, full), 1):
            rn = resname[(a['resn'], a['resi'])]
            rk = [93, 95, 118, 260, 261].index(a['resi']) + 1
            nm = a['name'] if len(a['name']) == 4 else ' ' + a['name']
            pdb.append(f"{'HETATM' if rk > 3 else 'ATOM  '}{k:5d} {nm:<4s} {rn:3s} A{rk:4d}    "
                       f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}  1.00  0.00          {a['el']:>2s}")
        pdb.append('END')
        (OUT / 'geom' / f'{pid}.pdb').write_text('\n'.join(pdb) + '\n')
        base = 'Rx' if (p['group'] == 'g' or pid.startswith('c2_') or pid == 'Rx') else 'X'
        sets = dict(gate_plan=pid in ['X', 'R'] or p['group'] in ('a', 'c'),
                    gate_rigid=pid in ['X', 'Rx'] or p['group'] in ('a', 'c'),
                    fit=pid not in ('R', 'N') and p['group'] in ('a', 'b', 'c', 'd'),
                    part_A=pid in ('X', 'Rx') or p['group'] in ('a', 'd', 'e', 'f', 'g'),
                    part_B=pid in ('X', 'Rx', 'Nx') or p['group'] == 'c',
                    report_only=pid in ('R', 'N'), base=base)
        listing.append(dict(order=n, id=pid, group=p['group'], definition=p['definition'],
                            parameters=p['parameters'], internal=p['internal'], contacts=p['contacts'],
                            flags=p['flags'], sets=sets,
                            xyz_sha256=hashlib.sha256((OUT / 'geom' / f'{pid}.xyz').read_bytes()).hexdigest(),
                            ligand_coordinates={nm: [round(float(t), 8) for t in v]
                                                for nm, v in zip(m.names, p['coords'])}))
    doc = dict(
        what='test set of gate G-R (md/PLAN.md, 6 Oct 2026, 3:02 AM IST, point 2); written by md/qm/gr_points.py',
        model=dict(file='md/qm/mcpb/CA2_small.pdb', natoms=64, fixed_atoms='1-46 (three ethyl-imidazoles, zinc)',
                   ligand_atoms='47-64', ligand_names=m.names, charge=1, multiplicity=1,
                   zinc=[float(t) for t in m.zn], thr199_og1=[float(t) for t in m.og1],
                   Zn_N2_crystal=m.d_zn_n, N2_S_crystal=m.d_n_s, N2_H=NH),
        limits=dict(model=LIMITS, protein=PROT_LIMITS, drop_heavy_heavy=DROP_HEAVY, grid_deg=GRID),
        like_points_check=like_check, path_a=path_info, omega_X=om_x, inversion_with_N2_fixed=n2_fixed,
        chosen_turns=chosen, turn_ranges_5deg_scan=ranges,
        added_points=dict(when='6 Oct 2026, after the review of the gate (md/PLAN.md, entry of 8:56 AM IST, point 2)',
                          e=e_ids, f=f_ids, g=g_ids, turns_about_Zn_N2_possible_multiples_of_15=zok,
                          turns_about_Zn_N2_chosen=zturns, omega_Rx=om_r),
        closest_cap_atom_over_all_points=min(p['internal']['closest_cap_atom'] for p in listing), dropped=[p for p in pts if 'drop' in pts[p]['flags']],
        n_points=len(listing), order=[p['id'] for p in listing], points=listing,
        versions=dict(numpy=np.__version__))
    text = json.dumps(doc, indent=1) + '\n'
    old = (OUT / 'points.json').read_text() if (OUT / 'points.json').exists() else None
    (OUT / 'points.json').write_text(text)

    cols = ['order', 'id', 'group', 'Zn_N2', 'N2_S', 'Zn_N2_S', 'Zn_N2_H', 'S_N2_H', 'angle_sum_N2', 'phi', 'omega',
            'Zn_N2_S_C4', 'H_N2_S_C4', 'H_N2_S_O1', 'H_N2_S_O2', 'C4_S_N2', 'N94_Zn_N2', 'N96_Zn_N2', 'N119_Zn_N2',
            'Zn_O1', 'Zn_O2', 'H_OG1', 'N2_H_OG1', 'rmsd_to_X_heavy', 'model_heavy_heavy', 'model_heavy_heavy_pair',
            'model_heavy_h', 'model_heavy_h_pair', 'model_h_h', 'model_h_h_pair', 'protein_polar_polar',
            'protein_polar_polar_pair', 'protein_other_heavy', 'protein_other_heavy_pair', 'protein_any_h',
            'protein_any_h_pair', 'flags', 'gate_plan', 'gate_rigid', 'fit', 'part_A', 'part_B', 'report_only', 'base',
            'N94_Zn_N2_S', 'signed_volume']
    rows = [','.join(cols)]
    for p in listing:
        c, po = p['contacts']['model'], p['contacts']['protein_outside_model']
        r = dict(p['internal'], order=p['order'], id=p['id'], group=p['group'],
                 model_heavy_heavy=c['heavy_heavy']['d'], model_heavy_heavy_pair=c['heavy_heavy']['pair'],
                 model_heavy_h=c['heavy_h']['d'], model_heavy_h_pair=c['heavy_h']['pair'],
                 model_h_h=c['h_h']['d'], model_h_h_pair=c['h_h']['pair'],
                 protein_polar_polar=po['polar_polar']['d'] if po['polar_polar'] else '',
                 protein_polar_polar_pair=po['polar_polar']['pair'] if po['polar_polar'] else '',
                 protein_other_heavy=po['other_heavy']['d'] if po['other_heavy'] else '',
                 protein_other_heavy_pair=po['other_heavy']['pair'] if po['other_heavy'] else '',
                 protein_any_h=po['any_h']['d'] if po['any_h'] else '',
                 protein_any_h_pair=po['any_h']['pair'] if po['any_h'] else '',
                 flags=' '.join(p['flags']), **p['sets'])
        rows.append(','.join(str(r[k]) for k in cols))
    (OUT / 'points.csv').write_text('\n'.join(rows) + '\n')
    sc = ['config,delta,Zn_N2_S_C4,H_N2_S_C4,heavy_heavy,pair_heavy_heavy,heavy_h,h_h,Zn_O1,Zn_O2,'
          'protein_polar_polar,protein_other_heavy,flags,chosen']
    for r in scan_rows:
        sc.append(','.join(str(t) for t in (
            r['config'], r['delta'], r['Zn_N2_S_C4'], r['H_N2_S_C4'], r['heavy_heavy'], r['pair_heavy_heavy'],
            r['heavy_h'], r['h_h'], r['Zn_O1'], r['Zn_O2'], r['protein_polar_polar'], r['protein_other_heavy'],
            ' '.join(r['flags']),
            int(r['delta'] in chosen[r['config']]))))
    (OUT / 'scan_ns_contacts.csv').write_text('\n'.join(sc) + '\n')
    zc = ['delta,N94_Zn_N2_S,heavy_heavy,pair_heavy_heavy,heavy_h,pair_heavy_h,h_h,pair_h_h,Zn_O1,Zn_O2,H_OG1,'
          'protein_polar_polar,protein_other_heavy,flags,chosen']
    for r in zscan:
        zc.append(','.join(str(t) for t in (
            r['delta'], r['N94_Zn_N2_S'], r['heavy_heavy'], r['pair_heavy_heavy'], r['heavy_h'], r['pair_heavy_h'],
            r['h_h'], r['pair_h_h'], r['Zn_O1'], r['Zn_O2'], r['H_OG1'], r['protein_polar_polar'],
            r['protein_other_heavy'], ' '.join(r['flags']), int(r['delta'] in zturns))))
    (OUT / 'scan_znn_contacts.csv').write_text('\n'.join(zc) + '\n')

    print(f"{len(listing)} points, order: {' '.join(doc['order'])}")
    print('turns chosen:', chosen, ' dropped:', doc['dropped'])
    for p in listing:
        i, c = p['internal'], p['contacts']
        po = c['protein_outside_model']['other_heavy']
        pp = c['protein_outside_model']['polar_polar']
        print(f"{p['order']:2d} {p['id']:8s} ZnN {i['Zn_N2']:.3f} ZnNS {i['Zn_N2_S']:6.1f} ZnNH {i['Zn_N2_H']:6.1f} "
              f"SNH {i['S_N2_H']:6.1f} sum {i['angle_sum_N2']:6.1f} phi {i['phi']:6.1f} om {i['omega']:6.1f} "
              f"tauZn {i['Zn_N2_S_C4']:7.1f} tauH {i['H_N2_S_C4']:7.1f} | model {c['model']['heavy_heavy']['d']:.2f} "
              f"{c['model']['heavy_h']['d']:.2f} {c['model']['h_h']['d']:.2f} | protein "
              f"{pp['d'] if pp else 9.99:.2f} {po['d'] if po else 9.99:.2f} {po['pair'] if po else ''} | H..OG1 {i['H_OG1']:.2f} {i['N2_H_OG1']:.0f} | {' '.join(p['flags'])}")
    if old is not None:
        print('points.json identical to the one of the previous run:', old == text)


if __name__ == '__main__':
    main()
