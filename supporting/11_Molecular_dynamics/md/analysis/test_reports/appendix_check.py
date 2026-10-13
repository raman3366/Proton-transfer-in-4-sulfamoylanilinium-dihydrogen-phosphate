# throwaway (version 1.3): is appendix A of SPEC.md enough to read the driver's files, and section 9D enough to make
# the count?  A small reader written from the TEXT of appendix A (DCD layout, PDB columns, manifest keys, time axis)
# and of sections 2, 9D and 9E, with the Python standard library and NumPy only: no MDAnalysis, nothing imported from
# analyse_md.py or run_md.py.  It makes the two water counts of section 9D for every full frame of a run
# (equilibration and production) and sets them, and the times, beside the CSV files that analyse_md.py wrote.
# usage: appendix_check.py OUT.json RUN_DIR ANALYSIS_DIR [RUN_DIR ANALYSIS_DIR ...]
import csv, json, math, os, re, struct, sys
import numpy as np


def read_dcd(path):
    """Appendix A.2: header of 276 bytes, then frames of 80 + 12 NATOM bytes; the number of frames from the file size."""
    size = os.path.getsize(path)
    with open(path, 'rb') as f:
        raw = f.read()
    if size < 276:
        return None
    assert struct.unpack('<i', raw[0:4])[0] == 84 and raw[4:8] == b'CORD'
    nset, istart, nsavc, _last = struct.unpack('<4i', raw[8:24])
    delta = struct.unpack('<f', raw[44:48])[0]
    has_cell = struct.unpack('<i', raw[48:52])[0]
    assert struct.unpack('<i', raw[88:92])[0] == 84 and struct.unpack('<i', raw[92:96])[0] == 164
    natom = struct.unpack('<i', raw[268:272])[0]
    assert has_cell == 1
    flen = 80 + 12 * natom
    n = (size - 276) // flen
    cells, xyz = np.zeros((n, 6)), np.zeros((n, natom, 3))
    for k in range(n):
        o = 276 + k * flen
        assert struct.unpack('<i', raw[o:o + 4])[0] == 48
        cells[k] = struct.unpack('<6d', raw[o + 4:o + 52])
        o += 56
        for d in range(3):
            assert struct.unpack('<i', raw[o:o + 4])[0] == 4 * natom
            xyz[k, :, d] = np.frombuffer(raw, dtype='<f4', count=natom, offset=o + 4)
            o += 8 + 4 * natom
    return {'nset': nset, 'istart': istart, 'nsavc': nsavc, 'timestep_ps': delta * 0.04888821, 'natom': natom, 'n_frames': n, 'cells': cells, 'xyz': xyz}


def lattice(cell):
    """Appendix A.2: the unit-cell record holds a, cos gamma, b, cos beta, cos alpha, c."""
    a, cg, b, cb, ca, c = cell
    sg = math.sqrt(1.0 - cg * cg)
    cx, cy = c * cb, c * (ca - cb * cg) / sg
    return np.array([[a, 0.0, 0.0], [b * cg, b * sg, 0.0], [cx, cy, math.sqrt(c * c - cx * cx - cy * cy)]])


def pdb_atoms(path):
    """Appendix A.1: every line that begins with ATOM or HETATM; atom name in columns 13 to 16, coordinates 31 to 54."""
    names, xyz = [], []
    for line in open(path):
        if line.startswith(('ATOM', 'HETATM')):
            names.append(line[12:16].strip())
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return names, np.array(xyz)


def fit(x, r):
    """Section 2, step 2: least squares, equal weights, rotation without reflection.  F(p) = (p - gx) Q + gr."""
    gx, gr = x.mean(0), r.mean(0)
    u, _, vt = np.linalg.svd((x - gx).T @ (r - gr))
    q = u @ np.diag([1.0, 1.0, np.sign(np.linalg.det(u @ vt))]) @ vt
    return q, gx, gr


def min_image(d, B):
    """Section 2: the smallest |d + i a + j b + k c|; fractional rounding first, then the 27 neighbouring cells."""
    d = d - np.rint(d @ np.linalg.inv(B)) @ B
    best = np.full(d.shape[:-1], np.inf)
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                best = np.minimum(best, np.linalg.norm(d + i * B[0] + j * B[1] + k * B[2], axis=-1))
    return best


def counts(frame, B, sel, ref_ca, ref_site, cut):
    q, gx, gr = fit(frame[sel['protein_CA']], ref_ca)
    site = (ref_site - gr) @ q.T + gx                       # section 9D, step 2: s_a = (r_a - g_r) Q^T + g_x
    w = frame[sel['water_O']]
    dw = np.min([min_image(w - s, B) for s in site], axis=0)
    return int((dw <= cut).sum()), float(dw.min())


def check_run(run_dir, ana_dir):
    man = json.load(open(f'{run_dir}/manifest.json'))
    sysdir = man['system_dir']
    sel = json.load(open(f'{sysdir}/selections.json'))
    rn, rx = pdb_atoms(f'{sysdir}/reference.pdb')
    ref_ca = rx[[i for i, n in enumerate(rn) if n == 'CA']]
    st = (man.get('hydration_check') or {}).get('settings')
    P, c, source = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5, 'specification (defaults)'
    if st is not None:                                       # section 9D, settings (the usable case only; these runs have it or lack the key)
        P, c, source = list(st['site_atoms']), float(st['water_cutoff_A']), 'manifest hydration_check.settings'
    site_p = np.array([rx[rn.index(a)] for a in P])
    site_w = np.array([rx[rn.index(a)] for a in ('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1')])
    n_full_pdb = len(pdb_atoms(f'{run_dir}/full.pdb')[0])
    dtf = man['protocol']['output']['full_frame_interval_ps']
    stages = [(s['tag'], s['time_start_ps'], s['time_end_ps']) for s in man['equilibration']['stages']]
    t_eq = max(e for _, _, e in stages)
    rows, head_ok, nset_ok = {'equil': [], 'prod': []}, True, True
    for tag, s0, s1 in stages:                               # the order of the manifest is the order of the run (appendix A.4)
        m = re.fullmatch(r'eq(\d+)(?:x([1-9]\d*))?', tag)
        d = read_dcd(f'{run_dir}/full_{tag}.dcd')
        if d is None or d['n_frames'] == 0:
            continue
        assert d['natom'] == n_full_pdb
        head_ok &= abs(d['istart'] * d['timestep_ps'] - (s0 + dtf)) <= 0.5 * d['timestep_ps'] + 1e-6 and abs(d['nsavc'] * d['timestep_ps'] - dtf) < 1e-6 * dtf + 1e-9
        nset_ok &= d['nset'] == d['n_frames']
        for k in range(d['n_frames']):
            B = lattice(d['cells'][k])
            npl, dpl = counts(d['xyz'][k], B, sel, ref_ca, site_p, c)
            nwd, dwd = counts(d['xyz'][k], B, sel, ref_ca, site_w, 3.5)
            rows['equil'].append((int(m.group(1)), int(m.group(2) or 0), (s0 + (k + 1) * dtf - t_eq) / 1000.0, npl, dpl, nwd, dwd))
    j = 0
    for seg in man['production']['segments']:
        d = read_dcd(f'{run_dir}/full_{seg["tag"]}.dcd')
        head_ok &= abs(d['istart'] * d['timestep_ps'] - (seg['time_start_ps'] + dtf)) <= 1e-6 and abs(d['nsavc'] * d['timestep_ps'] - dtf) < 1e-6 * dtf + 1e-9
        nset_ok &= d['nset'] == d['n_frames'] == seg['frames']['full']
        for k in range(d['n_frames']):
            B = lattice(d['cells'][k])
            npl, dpl = counts(d['xyz'][k], B, sel, ref_ca, site_p, c)
            nwd, dwd = counts(d['xyz'][k], B, sel, ref_ca, site_w, 3.5)
            rows['prod'].append((None, None, (j + 1) * dtf / 1000.0, npl, dpl, nwd, dwd))
            j += 1
    out = {'run_dir': run_dir, 'analysis_dir': ana_dir, 'site_atoms': P, 'water_cutoff_A': c, 'source': source,
           'n_atoms_full': n_full_pdb, 'n_water_O': len(sel['water_O']), 'n_CA': len(sel['protein_CA']),
           'dcd_headers_put_the_first_frame_one_interval_after_the_start_of_the_step': bool(head_ok),
           'nset_of_every_header_equals_the_frames_found_by_file_size': bool(nset_ok)}
    summ = json.load(open(f'{ana_dir}/summary.json'))
    wp = summ['water_ring_place']
    out['settings_as_the_analysis'] = wp['site_atoms'] == P and abs(wp['water_cutoff_A'] - c) < 1e-9 and wp['source'] == source
    worst, ok = 0.0, out['settings_as_the_analysis'] and head_ok and nset_ok
    for kind, name in (('equil', 'timeseries_equil_full.csv'), ('prod', 'timeseries_full.csv')):
        theirs = list(csv.DictReader(open(f'{ana_dir}/{name}')))
        mine = rows[kind]
        same_n = len(mine) == len(theirs)
        cnt = sum(1 for a, b in zip(mine, theirs) if a[3] == int(b['n_water_ring_place']) and a[5] == int(b['n_water_ring_wide_3p5A']))
        dd = max([max(abs(a[4] - float(b['d_water_ring_place_min_A'])), abs(a[6] - float(b['d_water_ring_wide_min_A']))) for a, b in zip(mine, theirs)] or [0.0])
        tt = max([abs(a[2] - float(b['time_ns'])) for a, b in zip(mine, theirs)] or [0.0])
        se = all((a[0], a[1]) == (int(b['stage']), int(b['extension'])) for a, b in zip(mine, theirs)) if kind == 'equil' else True
        worst = max(worst, dd)
        # the production series is written with three decimals of ns: a test run of sub-picosecond steps cannot be compared in time there
        t_ok = tt <= 1e-9 if kind == 'equil' else tt <= 0.0005 + 1e-12
        out[kind] = {'frames_here': len(mine), 'frames_of_the_analysis': len(theirs), 'frames_with_both_counts_identical': cnt,
                     'largest_difference_of_a_nearest_distance_A': round(dd, 6), 'largest_difference_of_a_time_ns': tt,
                     'stage_and_extension_identical': bool(se), 'counts_place_here': [a[3] for a in mine][:60]}
        ok &= same_n and cnt == len(mine) and dd <= 0.001 and t_ok and se
    out['largest_difference_of_a_nearest_distance_A'] = round(worst, 6)
    out['all_agree'] = bool(ok)
    return out


res = {'what': 'a reader written from appendix A and sections 2, 9D and 9E of SPEC.md alone (standard library and NumPy; no MDAnalysis, nothing '
               'of analyse_md.py or run_md.py): the two water counts of every full frame, and the times, against the files of analyse_md.py',
       'runs': [check_run(sys.argv[i], sys.argv[i + 1]) for i in range(2, len(sys.argv), 2)]}
res['n_full_frames'] = sum(r['equil']['frames_here'] + r['prod']['frames_here'] for r in res['runs'])
res['largest_difference_of_a_nearest_distance_A'] = max(r['largest_difference_of_a_nearest_distance_A'] for r in res['runs'])
res['all_agree'] = all(r['all_agree'] for r in res['runs'])
json.dump(res, open(sys.argv[1], 'w'), indent=1)
for r in res['runs']:
    print(os.path.basename(os.path.dirname(r['run_dir'])) + '/' + os.path.basename(r['run_dir']), '|', r['source'], r['site_atoms'], r['water_cutoff_A'], '| equilibration frames', r['equil']['frames_here'],
          'identical', r['equil']['frames_with_both_counts_identical'], '| production frames', r['prod']['frames_here'], 'identical', r['prod']['frames_with_both_counts_identical'],
          '| largest distance difference', r['largest_difference_of_a_nearest_distance_A'], 'A | headers', r['dcd_headers_put_the_first_frame_one_interval_after_the_start_of_the_step'],
          '| NSET', r['nset_of_every_header_equals_the_frames_found_by_file_size'], '| agree', r['all_agree'])
print('appendix check:', res['n_full_frames'], 'full frames,', 'ALL AGREE' if res['all_agree'] else 'DIFFERENCES')
