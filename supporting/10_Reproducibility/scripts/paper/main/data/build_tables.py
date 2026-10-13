#!/usr/bin/env python3
"""Build the main-text and SI tables (pandoc pipe tables with pandoc-crossref captions) from facts.json, the
derived/ JSON files and the result files.  Main: paper/main/tables/tableNN_*.md (labels as in the skeleton);
SI: paper/si/tables/tableSNN_*.md, numbered in SI-section order.  Include them in the manuscripts with
    <!-- include: tables/table02_proton.md -->
(paper/_shared/tools/include.lua).  Every printed number comes from facts.json or from the same files it cites.

    /root/workspace/cqc-dft/tools/paper-env/bin/python paper/main/data/build_tables.py
"""
import json
import re
import statistics
from decimal import Decimal, ROUND_HALF_EVEN
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RES = REPO / 'molecules/sfa/results'
DOCK = REPO / 'dock'
DER = HERE / 'derived'
MT = REPO / 'paper/main/tables'
ST = REPO / 'paper/si/tables'
MT.mkdir(parents=True, exist_ok=True)
ST.mkdir(parents=True, exist_ok=True)
HA_KCAL = 627.5095


def J(p):
    return json.loads(Path(p).read_text())


F = J(HERE / 'facts.json')
D1, D3, D4, D5, D6, D7, D8, D9 = (J(DER / f) for f in ['d01_xray.json', 'd03_geometry.json', 'd04_scan.json', 'd05_freq.json',
                                                       'd06_electronic.json', 'd07_hirshfeld.json', 'd08_periodic.json', 'd09_dock.json'])
D13 = J(DER / 'd13_referee_calcs.json')          # referee calculations (VM burst 3; facts sections 32-36)
INDEX = []
STORED_TIES = []          # (value, decimals) printed from a stored value that is itself a rounding tie (see f)


def f(x, d=3, sign=False, exact=False):
    """x at d decimals, rounded half to even on its full-precision value (the rounding rule recorded in facts.json
    _meta.rounding_rule; the digits of the float's shortest repr are taken as the value), true minus sign."""
    if x is None:
        return '–'
    if isinstance(x, str):
        return x
    tt = Decimal(repr(float(x))).as_tuple()
    if -tt.exponent == d + 1 and tt.digits[-1] == 5 and not exact:
        STORED_TIES.append((repr(x), d))           # a value already rounded to a tie: its full precision must be supplied
    q = Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_EVEN)
    s = f'{q:+.{d}f}' if sign else f'{q:.{d}f}'
    if s.startswith('-') and float(s) == 0:
        s = s[1:]
    return s.replace('-', '−')


def v(key):
    return F[key]['value']


def _rng(d_, n):
    """'min–max' at n decimals, or one value when both print the same (prose audit burst 3, Table S31)."""
    a_, b_ = f(d_['min'], n), f(d_['max'], n)
    return a_ if a_ == b_ else f'{a_}–{b_}'


def _auto_widths(header, align, rows):
    # 29 Sep 2026 (SI assembly): pandoc sizes the columns of a wide pipe table by the dash counts of the separator
    # line, so all-equal ':---' separators squeeze long first columns (run names, ligand names) into overlapping
    # text in the PDF. Where the separator was not sized by hand (every entry <= 4 characters), scale the dashes to
    # the longest cell of each column (3-30), keeping the alignment colons. Cell contents are unchanged.
    # 29 Sep 2026 (ACS integration): Tables 3 (table04), S11, S13, S14, S20 (tableS18), S22 (tableS18b) and S31
    # (tableS27) are sized by hand so that every column is at least as wide as its widest unbreakable token (bond
    # labels such as N1–H1A···O5, symmetry positions), measured in the PDF table font; tables of 10 or more columns
    # print at 8 pt in the PDF (crossref.lua, manuscript.latex).
    if any(len(a) > 4 for a in align):
        return align
    out = []
    for j, a in enumerate(align):
        n = max([len(str(header[j]))] + [len(str(r[j])) for r in rows if j < len(r)])
        n = max(3, min(30, n))
        out.append((':' if a.startswith(':') else '') + '-' * n + (':' if a.endswith(':') else ''))
    return out


_NO_END = re.compile(r'(?:\b(?:e\.g|i\.e|vs|cf|ca|approx|Fig|Figs|Eq|Eqs|Ref|Refs|No|al)|(?<![\w])[A-Z]|\*[A-Z])$')


def split_title(caption):
    """ACS Omega (Appendix 2): a table title is one phrase or sentence, and the details go in footnotes under the
    table. -> (title, rest): the first sentence of the caption (a full stop outside brackets and math, not after an
    abbreviation or an initial such as *E. coli*, followed by a space and a capital, a digit or markup) and the
    remaining sentences, which table() prints as the first note under the table. No text is changed or dropped."""
    depth, math = 0, False
    for i, ch in enumerate(caption):
        if ch == '$':
            math = not math
        elif ch in '([':
            depth += 1
        elif ch in ')]':
            depth -= 1
        elif ch == '.' and depth == 0 and not math and caption[i + 1:i + 2] == ' ':
            nxt = caption[i + 2:i + 3]
            if nxt and not nxt.islower() and not _NO_END.search(caption[:i]):
                return _split_colon(caption[:i + 1], caption[i + 2:].strip())
    return _split_colon(caption, '')


def _split_colon(title, rest):
    """A first sentence of more than 25 words of the form 'Subject: details.' keeps the subject as the title; the
    details open the note."""
    if len(title.split()) <= 25:
        return title, rest
    depth = 0
    for i, ch in enumerate(title):
        if ch in '([':
            depth += 1
        elif ch in ')]':
            depth -= 1
        elif ch == ':' and depth == 0 and title[i + 1:i + 2] == ' ':
            head, tail = title[:i].rstrip() + '.', title[i + 2:].strip()
            tail = tail[0].upper() + tail[1:] if tail[:1].islower() else tail
            return head, (tail + ' ' + rest).strip()
    return title, rest


def table(path, header, align, rows, caption, label, notes=None, doc='main', section=''):
    align = _auto_widths(header, align, rows)
    caption, rest = split_title(caption)
    notes = ([rest] if rest else []) + list(notes or [])
    L = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(align) + '|']
    for r in rows:
        L.append('| ' + ' | '.join(str(c) for c in r) + ' |')
    L += ['', f'Table: {caption} {{#tbl:{label}}}']
    if notes:
        L += [''] + [' '.join(notes)]
    Path(path).write_text('\n'.join(L) + '\n')
    INDEX.append((doc, Path(path).name, label, section, caption.split('.')[0]))
    if STORED_TIES:
        print(f'STORED TIE in {Path(path).name}: {STORED_TIES}')
        STORED_TIES.clear()


# ================================================================================================= MAIN TABLES
# Table 1 -- models
M = [('Ion pair (M1)', 'C₆H₉N₂O₂S⁺ · H₂PO₄⁻ (anion at *x*, *y*, *z* − 1)', 'M1', 'Primary model: proton position, geometry, spectra, wavefunction analyses'),
     ('Ion pair, def2-TZVP', 'as M1', 'M1-TZ', 'Basis-set check (IEF-PCM)'),
     ('Cation (M2)', 'C₆H₉N₂O₂S⁺', 'M2', 'Fragment: interaction energy, TD-DFT, descriptors; docking sensitivity case'),
     ('Dihydrogen phosphate (M3)', 'H₂PO₄⁻', 'M3', 'Fragment: interaction energy'),
     ('Sulfanilamide (M4)', 'C₆H₈N₂O₂S', 'M4', 'Neutral species at pH 7.4; primary docking ligand; proton-transfer energetics'),
     ('Sulfanilamide, def2-TZVP', 'as M4', 'M4-TZ', 'Basis-set check (IEF-PCM)'),
     ('Sulfonamidate (M5)', 'C₆H₇N₂O₂S⁻', 'M5', 'Docking species for the zinc site; ligand of the simulations'),
     ('Phosphoric acid (M6)', 'H₃PO₄', 'M6', 'Proton-transfer energetics'),
     ('Urea (M7)', 'CH₄N₂O', 'M7', 'NLO reference (gas phase)')]
rows = []
for name, comp, m, use in M:
    rows.append([name, comp, f(v(f'model.{m}.charge'), 0, sign=v(f'model.{m}.charge') != 0), v(f'model.{m}.n_atoms'), v(f'model.{m}.nao'), use])
rows.append(['Five-pair cluster (M8)', '5 × (C₆H₉N₂O₂S⁺ · H₂PO₄⁻), cut from the periodic cells', '0', v('qtaim.xclu.n_atoms'),
             f"{v('qtaim.xclu.nao')} (6-31+G(d,p))", 'QTAIM of the nine crystal hydrogen bonds, salt and proton on O5 (gas phase)'])
rows.append(['Periodic model (M9)', '2 × (C₆H₉N₂O₂S⁺ · H₂PO₄⁻), *Pc*, *Z* = 2', '0', v('model.M9.n_atoms'), 'plane waves (70/560 Ry)',
             'Proton position in the crystal: relaxations and the constrained N1–H1A profile (PBE-D3(BJ), PBE0-D3(BJ))'])
table(MT / 'table01_models.md',
      ['Model', 'Composition', 'Charge', 'Atoms', 'Basis functions', 'Use'],
      # column widths by content (ACS integration review, 30 Sep 2026): the Use column holds sentences, the numeric
      # columns at most a few characters
      [':--------------', ':----------------------', '------:', '-----:', '---------------:', ':-------------------------'], rows,
      'Models, their charges and sizes, and what each is used for. Basis functions are spherical 6-311++G(d,p) '
      'unless another basis is stated; all molecular models are closed-shell singlets cut from the crystal coordinates '
      '(urea from a force-field conformer).', 'models',
      notes=['IEF-PCM, integral-equation-formalism polarisable continuum model; TD-DFT, time-dependent density functional '
             'theory; QTAIM, quantum theory of atoms in molecules; NLO, nonlinear optical. '
             'The ion pair uses 443 of its 445 functions after removal of two near-linear dependencies. The cluster (M8) '
             f"holds {v('qtaim.xclu.n_pairs')} ion pairs that contain all nine hydrogen bonds of the crystal, eight of them at the central pair (N1–H1B···O6 joins two outer pairs), cut from "
             f"relaxed periodic cell (i) and from the cell with H1A held on O5 at {f(v('protocol.qe_start_OH'), 2)} Å; it is used for bond descriptors, and its energy only for its sign."], section='2.1')  # prose audit b3r2: O5 geometry named

# Table 2 -- proton metric
ENV = [('X-ray (deposited H1A)', 'xray'), ('Gas phase (*ε* = 1)', 'gas'), ('IEF-PCM (*ε* = 4)', 'eps4'),
       ('IEF-PCM (*ε* = 8)', 'eps8'), ('IEF-PCM (*ε* = 20)', 'eps20'), ('IEF-PCM (water), salt minimum', 'pcm'),
       ('IEF-PCM (water), neutral minimum^a^', 'scanO'), ('Periodic model, salt^b^', 'periodic')]
ENV += [(f"Periodic, N1–H1A held at {f(t_, 2)} Å^c^", f'periodic.{k_}') for t_, k_ in
        zip(v('periodic.profile.targets'), ['p110', 'p125', 'p140', 'p155', 'p170'])]
ENERGY = {'pcm': '0', 'periodic': '0',
          'scanO': f"{f(v('pt.pcm.dE'), 2)} (Δ*G*~298~ {f(v('pt.pcm.dG298'), 2)}) kcal mol⁻¹"}
for k_ in ['p110', 'p125', 'p140', 'p155', 'p170']:
    ENERGY[f'periodic.{k_}'] = (f"{f(v(f'periodic.profile.{k_}.PBE_kJ'), 1, sign=True)} "
                                f"({f(v(f'periodic.profile.{k_}.PBE0_kJ'), 1, sign=True)}) kJ mol⁻¹")
rows = []
for lab, k in ENV:
    # Numbers audit round 1 (29 Sep 2026): the X-ray distances to H1A are printed to 0.01 Å, the precision of the CIF
    # (N1–H1A 0.87(3), H1A···O5 1.89(3) Å). The third decimal depends on the coordinate set (the 3-decimal _atom_site H
    # coordinates give δ −1.023 Å, the full-precision .res coordinates −1.018 Å) and is not supported by the data.
    dh = 2 if k == 'xray' else 3
    rows.append([lab, f(v(f'proton.{k}.d_N1_H1A'), dh), f(v(f'proton.{k}.d_H1A_O5'), dh), f(v(f'proton.{k}.delta'), dh, sign=True),
                 f(v(f'proton.{k}.d_N1_O5')), f(v(f'proton.{k}.P1O5')),
                 'N1' if v(f'proton.{k}.side') == 'N' else 'O5', ENERGY.get(k, '–')])
assert v('periodic.profile.no_neutral_min') is True and v('pcm.oform.n_imag') == 0 and v('freq.n_imag') == 0
assert v('periodic.iir1.final_sides') == ['N', 'N'] and v('periodic.iir1.start_sides') == ['O', 'N']
table(MT / 'table02_proton.md',
      ['Environment', '*d*(N1–H1A) (Å)', '*d*(H1A···O5) (Å)', '*δ* (Å)', '*d*(N1···O5) (Å)', '*d*(P1–O5) (Å)', 'Proton nearer', 'Energy relative to the salt'],
      # first column widened for its header word in the Word file (9 pt table text; review of 5 Oct 2026)
      [':--------------', '----------------:', '------------------:', '--------:', '-----------------:', '---------------:', ':---------', '-----------------:'], rows,
      'Position of the proton in each environment. Proton metric *δ* = *d*(N1–H1A) − *d*(H1A···O5), the N1···O5 distance '
      'and the energy relative to the salt minimum of the same model, for the ion pair in each environment '
      '(B3LYP-D3(BJ)/6-311++G(d,p)) and for the periodic cell (PBE-D3(BJ), formula unit 1; both formula units are '
      'identical). The *ε* = 4–20 optimisations started from the water '
      'optimum, the gas-phase one from the X-ray structure.', 'proton',
      notes=[f"^a^ End point of the relaxed proton scan (chain B), a minimum without imaginary frequencies; the barrier from the salt side "
             f"is {f(v('scan.barrier'), 2)} kcal mol⁻¹. Neutral − salt: Δ*E* = {f(v('pt.pcm.dE'), 2)}, Δ(*E* + ZPE) = {f(v('pt.pcm.dEZPE'), 2)}, "
             f"Δ*H*~298~ = {f(v('pt.pcm.dH298'), 2)} and Δ*G*~298~ = {f(v('pt.pcm.dG298'), 2)} kcal mol⁻¹ (harmonic, Table S6). "
             f"^b^ Relaxed from the X-ray structure; starts with H1A on O5 in both formula units (*δ* = {f(v('proton.periodic.ii.delta'), 3, sign=True)} Å "
             f"after relaxation) or in one only (*δ* = {f(v('periodic.iir1.delta'), 3, sign=True)} Å) relax to the same salt. "
             '^c^ Constrained profile: N1–H1A of both formula units held, every other atom relaxed in the fixed X-ray cell. Energies '
             'per formula unit, PBE-D3(BJ) with the PBE0-D3(BJ) value (single points at the same geometries) in parentheses. No '
             'point is a minimum: the energy rises at every step and the constraint force pushes H1A back to N1 at every point '
             f"(d*E*/d*d* {f(v('periodic.profile.min_slope'), 0)}–{f(v('periodic.profile.max_slope'), 0)} kJ mol⁻¹ Å⁻¹). An earlier "
             f"configuration with O5–H1A held at {f(v('protocol.qe_start_OH'), 2)} Å (*d*(N1–H1A) {f(v('proton.periodic.iirfix.d_N1_H1A'), 3)} Å) lies at "
             f"{f(v('periodic.iirfix.dE_kJ'), 1, sign=True)} ({f(v('periodic.hybrid.e100.best_kJ'), 1, sign=True)}) kJ mol⁻¹, "
             f"{f(v('periodic.profile.iirfix_excess'), 1)} kJ mol⁻¹ above the PBE-D3(BJ) profile at that N1–H1A distance (Table S15). "
             'X-ray N–H lengths are short, so the X-ray *δ* is a lower bound.'], section='3.1')

# Table 3 -- X-ray vs DFT by bond class
XR = D3['models']['xray']['bond_lengths']
CL = [('C–C (ring)', ['C1-C2', 'C1-C6', 'C2-C3', 'C3-C4', 'C4-C5', 'C5-C6']), ('C1–N1, C4–S1', ['C1-N1', 'C4-S1']),
      ('S1–O1, S1–O2', ['O1-S1', 'O2-S1']), ('S1–N2', ['N2-S1']), ('P1–O3, P1–O4 (P–OH)', ['P1-O3', 'P1-O4']),
      ('P1–O5, P1–O6', ['P1-O5', 'P1-O6']),
      ('All heavy-atom bonds', ['C1-C2', 'C1-C6', 'C2-C3', 'C3-C4', 'C4-C5', 'C5-C6', 'C1-N1', 'C4-S1', 'O1-S1', 'O2-S1', 'N2-S1',
                                'P1-O3', 'P1-O4', 'P1-O5', 'P1-O6'])]
MODS = [('IEF-PCM (water)', 'pcm'), ('Gas phase', 'gas'), ('Periodic model', 'periodic_i')]
PRETTY = lambda b: {'N2-S1': 'S1–N2', 'O1-S1': 'S1–O1', 'O2-S1': 'S1–O2'}.get(b, b.replace('-', '–'))
rows = []
import math
for cname, bl in CL:
    for mlab, mk in MODS:
        bl_m = D3['models'][mk]['bond_lengths']
        dev = [(b, bl_m[b] - XR[b]) for b in bl]
        mae = sum(abs(x) for _, x in dev) / len(dev)
        rmsd = math.sqrt(sum(x * x for _, x in dev) / len(dev))
        wb, wx = max(dev, key=lambda t: abs(t[1]))
        rows.append([cname if mlab == MODS[0][0] else '', mlab, f(mae), f(rmsd), f'{f(wx, 3, sign=True)} ({PRETTY(wb)})'])
for mlab, mk in MODS:
    a = D3['models'][mk]['angle_stats_heavy']
    rows.append(['Heavy-atom angles' if mlab == MODS[0][0] else '', mlab, f(a['mae'], 1) + '°', f(a['rmsd'], 1) + '°',
                 f"{f(D3['models'][mk]['angle_devs'][a['worst']], 1, sign=True)}° ({a['worst'].replace('-', '–')})"])
table(MT / 'table03_geom.md', ['Bond class', 'Model', 'MAE (Å)', 'RMSD (Å)', 'Largest deviation (Å)'],
      [':-----------', ':------', '-------:', '--------:', ':------------------'], rows,
      'Deviations of the optimised bond lengths (and of the heavy-atom bond angles, last block) from the X-ray '
      'structure, by bond class. Ion pair at B3LYP-D3(BJ)/6-311++G(d,p); periodic model at PBE-D3(BJ) with the cell '
      'fixed at the X-ray values. Deviation = model − X-ray.', 'geom',
      notes=['Bonds to hydrogen are omitted: the X-ray C–H are riding (0.93 Å) and the O–H restrained (0.82 Å). '
             'The def2-TZVP and dispersion-free models are compared in Table S11.'], section='3.2 (printed as Table 4)')

# Table 4 -- hydrogen bonds
hb_x = {f"{h['D']}-{h['H']}...{h['A']}": h for h in D1['hbonds']}
hb_p = {h['name']: h for h in D8['structures']['relax_i']['hbonds']}
qp = {h['name']: h for h in D6['qtaim']['sfa_pair_pcm_d3']['hbonds']}
SYMTXT = {'1_554': '*x*, *y*, *z* − 1', '2_554': '*x*, −*y*, *z* − ½', '2_654': '*x* + 1, −*y*, *z* − ½', '2_564': '*x*, −*y* + 1, *z* − ½',
          '1_655': '*x* + 1, *y*, *z*', '1_455': '*x* − 1, *y*, *z*'}
pcm_geo = D3['models']['pcm']['proton']
rows = []
for key, h in hb_x.items():
    nm = key.replace('...', '···').replace('-', '–')
    per = hb_p[key]
    ip = qtx = None
    if key == 'N1-H1A...O5':
        ip, qtx = pcm_geo['d_N1_O5'], qp['N1-H1A...O5']
    elif key == 'C6-H6...O4':
        ip, qtx = pcm_geo['d_C6_O4'], qp['C6-H6...O4']
    kq = key.replace('-', '').replace('...', '_')
    rows.append([nm, SYMTXT[h['sym']], f(float(h['cif']['DA'].split('(')[0])), f(per['DA']), f(ip) if ip else '–',
                 f(float(h['cif']['DHA'].split('(')[0]), 0), f(per['DHA'], 1),
                 f(qtx['rho'], 4) if qtx else '–', f(qtx['lap'], 4) if qtx else '–',
                 f(v(f'qtaim.xclu.{kq}.salt.rho'), 4), f(v(f'qtaim.xclu.{kq}.d_rho_pct'), 0, sign=True) + ('^a^' if kq == 'N1H1A_O5' else '')])
table(MT / 'table04_hbonds.md',
      ['Contact', 'Acceptor at', '*d*(D···A) X-ray (Å)', '*d*(D···A) periodic (Å)', '*d*(D···A) ion pair (Å)',
       '∠D–H···A X-ray (°)', '∠D–H···A periodic (°)', '*ρ*, ion pair (a.u.)', '∇²*ρ*, ion pair (a.u.)',
       '*ρ*, cluster (a.u.)', 'Δ*ρ*, H1A on O5 (%)'],
      [':-----------------------------', ':-------------------------------', '----------------:', '----------------:', '----------------:', '----------------:', '----------------:', '-------------:', '-------------:', '-------------:', '------------:'], rows,
      'Hydrogen bonds of the crystal. Columns give the deposited geometry (X-ray), the periodic PBE-D3(BJ) model, the IEF-PCM (water) '
      'ion pair with its QTAIM descriptors at the H···A bond critical point for the two contacts it contains '
      '(B3LYP-D3(BJ)/6-311++G(d,p)), and *ρ* at the bond critical point of every contact in the five-pair cluster cut from '
      'the periodic salt, with its change when H1A sits on O5 '
      # prose audit burst 3 (re-review O3): name the O5 geometry, the frozen configuration, not a profile point
      f"(O5–H1A held at {f(v('protocol.qe_start_OH'), 2)} Å; B3LYP-D3(BJ)/6-31+G(d,p), gas phase).", 'hbonds',
      notes=['^a^ With H1A on O5 the critical point is that of N1···H1A. The cluster values are for the contact at the central '
             'ion pair, except N1–H1B···O6, which joins two outer pairs (Table S22); bond-energy estimates from the critical points are not given (Section S9). X-ray s.u.: '
             '0.002–0.003 Å on D···A, 3–5° on angles.'], section='3.1 (printed as Table 3)')

# Table 5 -- FMO and descriptors
FS = [('Ion pair', 'IEF-PCM (water)', 'pcm', 'pcm'), ('Ion pair (proton on O5)', 'Gas phase', 'gas', None),
      ('Cation', 'IEF-PCM (water)', 'cation.pcm', 'cation.pcm'), ('Cation', 'Gas phase', 'cation.gas', 'cation.gas'),
      ('Sulfanilamide', 'IEF-PCM (water)', 'sfa.pcm', 'sfa.pcm'), ('Sulfanilamide', 'Gas phase', 'sfa.gas', 'sfa.gas'),
      ('Sulfonamidate', 'IEF-PCM (water)', 'sfaan.pcm', None)]
rows = []
for sp, mod, fk, dk in FS:
    r_ = [sp, mod, f(v(f'fmo.{fk}.homo')), f(v(f'fmo.{fk}.lumo')), f(v(f'fmo.{fk}.gap'))]
    if dk:
        r_ += [f(v(f'dscf.{dk}.I')), f(v(f'dscf.{dk}.A')), f(v(f'dscf.{dk}.eta')), f(v(f'dscf.{dk}.chi')), f(v(f'dscf.{dk}.omega'))]
    else:
        r_ += ['–'] * 5
    rows.append(r_)
table(MT / 'table05_fmo.md',
      ['Species', 'Model', '*E*~HOMO~', '*E*~LUMO~', 'Δ*E*~gap~', '*I*', '*A*', '*η*', '*χ*', '*ω*'],
      [':--------', ':------', '---:', '---:', '---:', '---:', '---:', '---:', '---:', '---:'], rows,
      'Frontier orbital energies and vertical ΔSCF reactivity descriptors (eV), B3LYP-D3(BJ)/6-311++G(d,p). '
      '*I* = *E*(*N*−1) − *E*(*N*), *A* = *E*(*N*) − *E*(*N*+1), *η* = (*I* − *A*)/2, *χ* = (*I* + *A*)/2, *ω* = *χ*²/(2*η*).', 'fmo',
      notes=[f"Orbital-based (Koopmans-type) values are compared in Table S16. The gas-phase *A* of sulfanilamide is negative "
             f"({f(v('dscf.sfa.gas.A'), 3)} eV): the extra electron is unbound and that ΔSCF anion is a basis-set artefact. The gas-phase "
             f"cation is the anilinium tautomer, which is not the gas-phase minimum of the isolated cation [@uhlemann2021preferred]."], section='3.5')

# Table 6 -- NLO
NL = [('Ion pair', 'IEF-PCM (water)', 'pcm'), ('Ion pair (IEF-PCM geometry)', 'Gas phase', 'gasAtPcm'),
      ('Sulfanilamide', 'IEF-PCM (water)', 'sfa'), ('Urea', 'Gas phase', 'urea')]
rows = [[sp, mod, f(v(f'nlo.{k}.mu'), 2), f(v(f'nlo.{k}.alpha'), 1), f(v(f'nlo.{k}.beta'), 1), f(v(f'nlo.{k}.ratio_urea'), 1)] for sp, mod, k in NL]
table(MT / 'table06_nlo.md', ['Species', 'Model', '*μ* (D)', '*ᾱ* (a.u.)', '*β*~tot~ (a.u.)', '*β*~tot~/*β*~tot~(urea)'],
      [':--------', ':------', '---:', '---:', '---:', '---:'], rows,
      'Static dipole moment, isotropic polarisability and first hyperpolarisability, B3LYP-D3(BJ)/6-311++G(d,p), finite '
      'field (0.001 and 0.002 a.u., Richardson-extrapolated), with urea at the same level. *β* in the Taylor (T) '
      'convention; 1 a.u. = 8.639 × 10⁻³³ esu; the B convention halves every *β*.', 'nlo',
      notes=['Electronic, static values of one ion pair: no vibrational or frequency-dependent terms, no local-field '
             'factors. The like-for-like comparison with urea is the gas-phase row.'], section='not printed since the revision (values in SI Table S27)')

# Table 7 -- redocking / cross-docking
vs = J(DOCK / 'analysis/summary.json')
RD = [('7S2K (Sul2)', 'Dihydropteroate (redock)', 'Vina', 'dock.val.T5_7S2K.dhp_anion__vina', 'validation'),
      ('7S2K (Sul2)', 'Dihydropteroate (redock)', 'Vinardo', 'dock.val.T5_7S2K.dhp_anion__vinardo', 'validation'),
      ('3HS4 (hCA II)', 'Sulfanilamide (cross-dock)', 'AD4Zn', 'dock.val.T4_3HS4.sfa_neutral__ad4zn', 'validation'),
      ('6RL9 (hCA II)', 'Sulfanilamide (redock)', 'AD4Zn', 'dock.val.T3_6RL9.sfa_neutral__ad4zn', 'validation'),
      ('6RL9 (hCA II)', 'Sulfonamidate (redock)', 'AD4Zn', 'dock.val.T3_6RL9.sfa_anion__ad4zn', 'validation'),
      ('6RL9 (hCA II)', 'Sulfanilamide (redock)', 'Vina', 'dock.val.T3_6RL9.sfa_neutral__vina', 'validation'),
      ('3TZF (DHPS, chain B)', 'Sulfamethoxazole anion (redock)', 'Vina', 'dock.val.T1_3TZF.smx_anion__vina', 'validation'),
      ('3HS4 (hCA II)', 'Acetazolamide (redock)', 'AD4Zn', 'dock.val.T4_3HS4.azm_neutral__ad4zn', 'validation'),
      ('1AJ0 (DHPS, *E. coli*)', 'Sulfanilamide (redock)', 'Vina', 'dock.val.T2_1AJ0.sfa_neutral__vina', 'validation')]
rows = []
for rec, lig, prog, k, stage in RD:
    rows.append([rec, lig, prog, 'ETKDG', f"{f(v(k + '.rmsd'), 2)} ± {f(v(k + '.rmsd_sd'), 2)}", f"{f(v(k + '.f1'), 1)} ({f(v(k + '.f2'), 1)})",
                 (v(k + '.G2').lower() + ' / ' + v(k + '.crit2A')) if k + '.G2' in F else '–'])
ps = J(DOCK / 'results/production_summary.json')
for rec, lig, prog, key_, job in [('6RL9 (hCA II)', 'Sulfanilamide (redock)', 'AD4Zn', 'dock.prod.6RL9.sfa.ad4zn.rmsd', 'T3_6RL9/sfa_neutral__ad4zn__prod'),
                                  ('3HS4 (hCA II)', 'Sulfanilamide (cross-dock)', 'AD4Zn', 'dock.prod.3HS4.sfa.ad4zn.rmsd', 'T4_3HS4/sfa_neutral__ad4zn__prod')]:
    j = ps['jobs'][job]
    full = [x for kk, x in j['rmsd_to_crystal'].items() if kk.startswith('SAN')][0]['full']
    assert j['posebusters']['pass_fraction'] == 1.0
    rows.append([rec, lig, prog, 'DFT (20 seeds)', f"{f(v(key_), 2)} ± {f(v(key_ + '_sd'), 2)}", f"{f(full['f_le_1A'], 1)} ({f(full['f_le_2A'], 1)})",
                 '–'])
assert v('dock.glu106.6RL9.sfa_anion.ad4zn.val.carb.gate') == 'pass' and v('dock.glu106.6RL9.sfa_anion.ad4zn.val.acid.gate') == 'fail'
table(MT / 'table07_redock.md',
      ['Receptor', 'Ligand', 'Scoring function', 'Ligand geometry', 'RMSD (Å)', 'Success fraction ≤ 1 Å (≤ 2 Å)', 'Strict gate / 2.0 Å criterion'],
      [':---', ':---', ':---', ':---', '---:', '---:', ':---'], rows,
      'Redocking and cross-docking of the crystal ligands. RMSD is the symmetry-corrected heavy-atom RMSD of the top-ranked '
      'pose in the receptor frame (mean ± SD over seeds, 10 seeds unless stated); the success fraction is the fraction of '
      'seeds whose top pose lies within 1 Å (2 Å). Every run passed the PoseBusters physics checks (pass fraction 1.0). Exhaustiveness 32, 20 Å box; '
      'AD4Zn, AutoDock4Zn.', 'redock',
      # 11 Oct 2026 (stage S7, row 29 of map_text.md): the note names the state of Glu106 in the hCA II receptors and
      # what the repeat with the carboxylate changes (facts dock.glu106.*; Table S36 = tableS31_dock_glu106_water).
      notes=['Validation used ETKDG/MMFF94s conformers; the last two rows repeat the hCA II sulfanilamide runs with the DFT '
             'geometry used in production. The hCA II rows are for the receptors as prepared, with Glu106 as the neutral '
             'acid (Section 2.5); with Glu106 as a carboxylate the 6RL9 run from the DFT geometry gives '
             f"{f(v('dock.glu106.6RL9.sfa_neutral.ad4zn.prod.carb.rmsd'), 2)} ± {f(v('dock.glu106.6RL9.sfa_neutral.ad4zn.prod.carb.rmsd_sd'), 2)} Å "
             f"with {v('dock.glu106.6RL9.sfa_neutral.ad4zn.prod.carb.n_le_1')} of {v('dock.glu106.6RL9.sfa_neutral.ad4zn.prod.carb.n_seeds')} "
             'seeds within 1 Å, and the sulfonamidate redock then passes the strict gate (Table S36). '
             'The 6RL9 AD4Zn result is bimodal in validation (near-native pose found in every '
             f"seed but ranked first in half); in that run the two orientations differ by {f(v('dock.val.6RL9.ad4zn.orientation_sep')[0], 3)}–"
             f"{f(v('dock.val.6RL9.ad4zn.orientation_sep')[1], 3)} kcal mol⁻¹ (over eight AD4Zn runs the near-native pose trails the top pose by "
             f"{f(v('dock.ad4zn8.orientation_gap')[0], 3)}–{f(v('dock.ad4zn8.orientation_gap')[1], 3)}; near-native: within 1.0 Å of the crystal pose), "
             'and which one ranks first depends on the input geometry and on the state of Glu106. The Vina run in 6RL9 misses the strict gate by its RMSD '
             'and fails the 2.0 Å criterion on the zinc contact '
             f"(Zn–N {f(v('dock.val.T3_6RL9.sfa_neutral__vina.ZnN'), 2)} Å against the {f(v('dock.gate.ZnN_max'), 1)} Å of the site check), not on the distance from the crystal pose. "
             'The 3HS4 cross-dock reference is 6RL9 sulfanilamide carried over by Cα '
             'superposition (0.328 Å).'], section='3.6')

# Table 8 -- scores and ranks
pc = J(DOCK / 'results/production_contacts.json')
SPN = {'sfa_neutral': 'Sulfanilamide', 'sfa_anion': 'Sulfonamidate', 'sfa_cation': '4-Sulfamoylanilinium', 'smx_neutral': 'Sulfamethoxazole',
       'smx_anion': 'Sulfamethoxazole anion', 'azm_neutral': 'Acetazolamide', 'azm_anion': 'Acetazolamide anion',
       'paba_anion': '*p*-Aminobenzoate', 'dhp_anion': '7,8-Dihydropteroate (crystal ligand, 78H)'}
REC = [('7S2K (Sul2)', 'T5_7S2K', 'vina', 'Vina'), ('3TZF (DHPS)', 'T1_3TZF', 'vina', 'Vina'),
       ('6RL9 (hCA II)', 'T3_6RL9', 'ad4zn', 'AD4Zn'), ('3HS4 (hCA II)', 'T4_3HS4', 'ad4zn', 'AD4Zn')]
rows = []
for rlab, t, sc, sclab in REC:
    ro = ps['rank_order'][f'{t}/{sc}']['ranking']
    first = True
    for rr in ro:
        job = f"{t}/{rr['species']}__{sc}__prod"
        c1 = pc.get(job, {}).get('1', {})
        def rname(k):
            if k.startswith('ZN'):
                return 'Zn'
            if k.startswith('HH2'):
                return 'pterin (HH2)'
            if k.startswith('POP'):
                return 'PPi (POP)'
            if k.startswith('HOH'):
                return 'water'
            if k.startswith('MG'):
                return 'Mg'
            return k[:1] + k[1:3].lower() + k[3:]
        res = ', '.join(rname(k) for k in list(c1.get('residues_within_4A', {}).keys())[:4])
        name = SPN[rr['species']] + (' (neutral form of the crystal ligand 08D)' if (t == 'T1_3TZF' and rr['species'] == 'smx_neutral') else '')
        bk = f"bulk.dockprod.{t}.{sc}.{rr['species']}"
        # ties: species whose mean scores lie within 0.1 kcal/mol share a range of ranks (referee chem m3, dock minor 1)
        means = [v(f"bulk.dockprod.{t}.{sc}.{x['species']}.score_mean") for x in ro]
        me = v(bk + '.score_mean')
        lo = 1 + sum(1 for x in means if x < me - 0.1 + 1e-9)
        hi = sum(1 for x in means if x < me + 0.1 - 1e-9)
        rk = f'{lo}/{len(ro)}' if lo == hi else f'{lo}–{hi}/{len(ro)}'
        rows.append([f'{rlab}, {sclab}' if first else '', name, f"{f(me, 1, exact=True)} ± {f(v(bk + '.score_sd'), 2)}",
                     rk, res])
        if rr['species'] == 'sfa_neutral':
            tk = f"dock.ties.{t.split('_')[1]}.{sc}"
            if tk in F:
                assert [lo, hi] == list(F[tk]['value']), (tk, lo, hi, F[tk]['value'])
        first = False


def _glu106_carb_score():
    """Mean AD4Zn score of sulfanilamide (DFT geometry, 20 seeds) with Glu106 as a carboxylate, at the one decimal of
    this table; the note prints one number for both hCA II structures, so the two must print alike, and acetazolamide
    must score better than sulfanilamide in both validation runs (stage S7, 11 Oct 2026)."""
    a = f(v('dock.glu106.6RL9.sfa_neutral.ad4zn.prod.carb.score'), 1, exact=True)
    b = f(v('dock.glu106.3HS4.sfa_neutral.ad4zn.prod.carb.score'), 1, exact=True)
    assert a == b, (a, b)
    assert v('dock.glu106.6RL9.azm_minus_sfa.carb') > 0 and v('dock.glu106.3HS4.azm_minus_sfa.carb') > 0
    return a


table(MT / 'table08_scores.md', ['Receptor, scoring', 'Ligand', 'Score (kcal mol⁻¹)', 'Rank', 'Nearest residues (≤ 4 Å)'],
      [':---', ':---', '---:', '---:', ':---'], rows,
      'Production docking scores and ranks. Scores are means ± SD over 10 seeds (20 for sulfanilamide with AD4Zn); ranks '
      'are by score among the species docked to each receptor, and the residues are those closest to the top pose of the '
      'largest cluster. The sulfanilamide '
      'species were docked from their IEF-PCM DFT geometries, the controls from ETKDG/MMFF94s conformers. Scores are '
      'empirical docking scores, not free energies of binding; only their order is interpreted. AD4Zn, AutoDock4Zn.', 'scores',
      notes=['Species whose mean scores lie within 0.1 kcal mol⁻¹ of each other share a range of ranks; the SD measures the '
             'reproducibility of the search over seeds, not the error of the scoring function. Ranks carry no site filter: '
             'with AD4Zn, sulfamethoxazole ranks first at both hCA II structures although its top pose does not reach the '
             'zinc (Zn–N 3.28–3.56 Å). In 3TZF the pose-validated reference is the sulfamethoxazole anion. PPi = pyrophosphate, '
             'HH2 = 6-hydroxymethylpterin diphosphate (pterin pyrophosphate). Scores of 10–23-heavy-atom ligands are size-biased '
             '(Section S13). '
             # 11 Oct 2026 (stage S7, row 30 of map_text.md): the state of Glu106 behind the hCA II rows
             'The hCA II rows are for the receptors with Glu106 as the neutral acid; with the carboxylate sulfanilamide scores '
             f"{_glu106_carb_score()} kcal mol⁻¹ in both structures, and acetazolamide stays ahead of it in the validation runs "
             '(Tables S36 and S37).'], section='3.7')

# ================================================================================================= SI TABLES
# S1 settings
CASES = [('sfa_pair_pcm_d3', 'Ion pair', 'B3LYP-D3(BJ)', '6-311++G(d,p)', 'IEF-PCM (water)'),
         ('sfa_pair_gas_d3', 'Ion pair', 'B3LYP-D3(BJ)', '6-311++G(d,p)', 'Gas phase'),
         ('sfa_pair_pcm_b3lyp', 'Ion pair', 'B3LYP', '6-311++G(d,p)', 'IEF-PCM (water)'),
         ('sfa_pair_pcm_d3_tz', 'Ion pair', 'B3LYP-D3(BJ)', 'def2-TZVP', 'IEF-PCM (water)'),
         ('sfa_pair_pcm_d3 (eps_ladder)', 'Ion pair', 'B3LYP-D3(BJ)', '6-311++G(d,p)', 'IEF-PCM, ε = 4, 8, 20'),
         ('sfa_pair_pcm_d3 (scan)', 'Ion pair', 'B3LYP-D3(BJ)', '6-311++G(d,p)', 'IEF-PCM (water), relaxed scan'),
         ('sfa_pair_pcm_d3 (fsp)', 'Ion pair', 'ωB97X-D, M06-2X, PBE0-D3(BJ), DF-MP2', '6-311++G(d,p)', 'IEF-PCM (water), single points'),
         ('sfa_pair_pcm_d3 (TD)', 'Ion pair', 'CAM-B3LYP, B3LYP', '6-311++G(d,p)', 'IEF-PCM (water), non-equilibrium'),
         ('cation / h2po4 / h3po4 / sulfanilamide _pcm_d3, _gas_d3', 'Fragments', 'B3LYP-D3(BJ)', '6-311++G(d,p)', 'IEF-PCM (water) and gas'),
         ('sulfonamidate_pcm_d3', 'Sulfonamidate', 'B3LYP-D3(BJ)', '6-311++G(d,p)', 'IEF-PCM (water)'),
         ('sulfanilamide_pcm_d3_tz', 'Sulfanilamide', 'B3LYP-D3(BJ)', 'def2-TZVP', 'IEF-PCM (water)'),
         ('urea_gas_d3', 'Urea', 'B3LYP-D3(BJ)', '6-311++G(d,p)', 'Gas phase'),
         ('periodic_qe', 'Periodic cell', 'PBE-D3(BJ)', 'SSSP 1.3.0, 70/560 Ry', 'Crystal, 2×1×1 k-points'),
         ('periodic_qe (hybrid)', 'Periodic cell', 'PBE0-D3(BJ), PBE-D3(BJ)', 'SG15 ONCV, 100/400 Ry (also 80/320)', 'Crystal, single points at the PBE-D3(BJ) geometries'),
         ('periodic_qe (profile)', 'Periodic cell', 'PBE-D3(BJ); PBE0-D3(BJ) single points', 'SSSP 1.3.0, 70/560 Ry; SG15 ONCV, 100/400 Ry', 'Crystal, N1–H1A held (damped dynamics)'),
         ('sfa_pair_pcm_d3_oform', 'Ion pair, neutral minimum', 'B3LYP-D3(BJ)', '6-311++G(d,p)', 'IEF-PCM (water), Hessian'),
         ('sfa_xclu_salt, sfa_xclu_neutral', 'Five-pair cluster (M8)', 'B3LYP-D3(BJ)', '6-31+G(d,p)', 'Gas phase, single points at periodic geometries'),
         ('sulfanilamide / sulfonamidate _pcm_d3 (torsion_scan)', 'Sulfanilamide, sulfonamidate', 'B3LYP-D3(BJ)', '6-311++G(d,p)', 'IEF-PCM (water), relaxed torsion scan')]
# SCF thresholds as run (dft/optimize.py, dft/proton_tools.py, dft/fsp.py, dft/tddft.py, dft/crystal_cluster_qtaim.py):
# optimisations and scans 1e-8 Ha in IEF-PCM and 1e-9 Ha in the gas phase; functional check and TD-DFT 1e-9 Ha;
# cluster single points 1e-8 Ha. (Until 29 Sep 2026 this column printed 1e-9 Ha for every molecular run.)
SCF_TH = {'sfa_pair_pcm_d3': '10⁻⁸ Ha', 'sfa_pair_gas_d3': '10⁻⁹ Ha', 'sfa_pair_pcm_b3lyp': '10⁻⁸ Ha', 'sfa_pair_pcm_d3_tz': '10⁻⁸ Ha',
          'sfa_pair_pcm_d3 (eps_ladder)': '10⁻⁸ Ha', 'sfa_pair_pcm_d3 (scan)': '10⁻⁸ Ha', 'sfa_pair_pcm_d3 (fsp)': '10⁻⁹ Ha',
          'sfa_pair_pcm_d3 (TD)': '10⁻⁹ Ha', 'cation / h2po4 / h3po4 / sulfanilamide _pcm_d3, _gas_d3': '10⁻⁸ Ha (IEF-PCM), 10⁻⁹ Ha (gas)',
          'sulfonamidate_pcm_d3': '10⁻⁸ Ha', 'sulfanilamide_pcm_d3_tz': '10⁻⁸ Ha', 'urea_gas_d3': '10⁻⁹ Ha', 'sfa_pair_pcm_d3_oform': '10⁻⁸ Ha',
          'sfa_xclu_salt, sfa_xclu_neutral': '10⁻⁸ Ha', 'sulfanilamide / sulfonamidate _pcm_d3 (torsion_scan)': '10⁻⁸ Ha'}


def _conv(c):
    return SCF_TH[c] if 'periodic' not in c else ('10⁻⁶ Ry' if 'hybrid' in c else '10⁻⁸ Ry' if 'profile' not in c else '10⁻⁹ Ry')


def _src(c):
    # Folders of the data package (06_DFT_results/, named in the caption; the same names as molecules/sfa/results/ of the
    # repository, periodic_qe/ packaged as periodic_QE/), and no <...>
    # placeholders: pandoc drops <species> as raw HTML and the | of <pcm|gas> split the table cell.
    if 'torsion_scan' in c:
        return '(species)_pcm_d3/ torsion_scan/'
    if 'profile' in c:
        return 'periodic_QE/ profile.json'
    if 'xclu' in c:
        return 'sfa_xclu_salt/, sfa_xclu_neutral/, crystal_cluster_qtaim/'
    if '/' in c:
        return '(fragment)_pcm_d3/, (fragment)_gas_d3/'
    return f'{c.split(" ")[0]}/'


rows = [[c, m, fx, b, env, ('level 3' if 'xclu' in c else 'level 4') if 'periodic' not in c else '–', _conv(c), _src(c)] for c, m, fx, b, env in CASES]
table(ST / 'tableS01_settings.md', ['Case', 'Model', 'Functional', 'Basis', 'Environment', 'Grid', 'SCF convergence', 'Source'],
      [':----------------', ':--------', ':-----------', ':-------------', ':-----------', ':-----', ':-------', ':-------------'], rows,
      'Settings of every calculation, with the case name and the folder of the data package (06_DFT_results/) that holds its results. PySCF 2.14.0 with density '
      'fitting and geomeTRIC 1.1.1 for the molecular models; IEF-PCM with *ε* = 78.3553 and Lebedev order 29; Quantum '
      'ESPRESSO 7.5 for the periodic model; Multiwfn 2026.9.20 for the wavefunction analyses.', 's-settings',
      notes=['The IEF-PCM ion pair and the dispersion-free pair were accepted on negligible forces (rms 7.35 × 10⁻⁵ and '
             '5.92 × 10⁻⁵ Ha Bohr⁻¹) after the job was stopped; the IEF-PCM pair is a confirmed minimum (0 imaginary '
             'frequencies). DF-MP2 is frozen-core with the aug-cc-pVTZ-RI fitting basis; its IEF-PCM term is the reaction field of the '
             'RHF density. The periodic PBE0-D3(BJ) runs use exact exchange with the adaptively compressed operator and a 1×1×1 '
             '(also 2×1×1) *q* mesh. Final single points for orbitals, Hessians, NLO and interaction energies used 10⁻¹⁰ Ha. The '
             'constrained periodic profile converged on the constrained forces below 10⁻³ Ry Bohr⁻¹ (damped dynamics; the '
             'unconstrained relaxations used BFGS and 5 × 10⁻⁴ Ry Bohr⁻¹).'], doc='si', section='S1')

# S2 energies of the optimised structures
OPT = ['sfa_pair_pcm_d3', 'sfa_pair_gas_d3', 'sfa_pair_pcm_b3lyp', 'sfa_pair_pcm_d3_tz', 'cation_pcm_d3', 'cation_gas_d3',
       'h2po4_pcm_d3', 'h2po4_gas_d3', 'h3po4_pcm_d3', 'h3po4_gas_d3', 'sulfanilamide_pcm_d3', 'sulfanilamide_gas_d3',
       'sulfanilamide_pcm_d3_tz', 'sulfonamidate_pcm_d3', 'urea_gas_d3']
# Names of the structures by model (Table 1) and environment; the result directory is kept as provenance (closing items,
# 29 Sep 2026: the SI named the structures by directory only).
STRUCT = {'sfa_pair_pcm_d3': 'Ion pair (M1), salt minimum, IEF-PCM (water)',
          'sfa_pair_pcm_d3_oform': 'Ion pair (M1), neutral (O5-protonated) minimum, IEF-PCM (water)',
          'sfa_pair_gas_d3': 'Ion pair (M1) optimised in the gas phase: neutral adduct',
          'sfa_pair_pcm_b3lyp': 'Ion pair (M1), IEF-PCM (water), B3LYP without D3(BJ)',
          'sfa_pair_pcm_d3_tz': 'Ion pair (M1), IEF-PCM (water), def2-TZVP',
          'cation_pcm_d3': 'Cation (M2), IEF-PCM (water)', 'cation_gas_d3': 'Cation (M2), gas phase',
          'h2po4_pcm_d3': 'Dihydrogen phosphate (M3), IEF-PCM (water)', 'h2po4_gas_d3': 'Dihydrogen phosphate (M3), gas phase',
          'h3po4_pcm_d3': 'Phosphoric acid (M6), IEF-PCM (water)', 'h3po4_gas_d3': 'Phosphoric acid (M6), gas phase',
          'sulfanilamide_pcm_d3': 'Sulfanilamide (M4), IEF-PCM (water)', 'sulfanilamide_gas_d3': 'Sulfanilamide (M4), gas phase',
          'sulfanilamide_pcm_d3_tz': 'Sulfanilamide (M4), IEF-PCM (water), def2-TZVP',
          'sulfonamidate_pcm_d3': 'Sulfonamidate (M5), IEF-PCM (water)', 'urea_gas_d3': 'Urea (M7), gas phase'}
rows = []
for c in OPT:
    o = J(RES / c / 'opt_summary.json')
    fr = RES / c / 'freq.json'
    nim = J(fr)['n_imaginary'] if fr.exists() else '–'
    zpe = f(J(fr)['ZPE_Ha'], 6) if fr.exists() else '–'
    rows.append([STRUCT[c], f(o['E_final_Ha'], 6), o['steps'], nim, zpe, f'`{c}`'])
    if c == 'sfa_pair_pcm_d3':
        # the neutral minimum: the chain-B anchor of the relaxed scan, with its own Hessian (VM burst 3)
        fo = J(RES / 'sfa_pair_pcm_d3_oform/freq.json')
        assert v('pcm.oform.same_as_scan_B') is True
        rows.append([STRUCT['sfa_pair_pcm_d3_oform'] + '^a^', f(fo['E_Ha'], 6), [p_['steps'] for p_ in D4['points'] if p_['chain'] == 'B' and p_['index'] == 0][0], fo['n_imaginary'], f(fo['ZPE_Ha'], 6), '`sfa_pair_pcm_d3_oform`'])
table(ST / 'tableS02_energies.md', ['Structure', '*E* (Ha)', 'Steps', 'Imaginary frequencies', 'ZPE (Ha)', 'Result directory'],
      [':-----', '---:', '---:', '---:', '---:', ':-----'], rows,
      'Energies of the optimised structures. Electronic energies (including the D3(BJ) term where used), the number of '
      'optimisation steps, and, where a Hessian was computed, the number of imaginary frequencies and the unscaled '
      'zero-point energy. Cartesian coordinates follow in Section S2.', 's-energies',
      notes=['"–": no Hessian was computed; those structures are stationary points on the forces criterion only. ^a^ Neutral '
             '(O5-protonated) minimum of the IEF-PCM ion pair: the end point of scan chain B (Section S5), whose Hessian was '  # prose audit burst 3: no forward table reference in a note
             'computed separately.'], doc='si', section='S2')
# coordinates (plain blocks)
L = ['<!-- Section S2 coordinates: generated by paper/main/data/build_tables.py from molecules/sfa/results/*/opt.xyz -->', '']
for c in OPT[:1] + ['sfa_pair_pcm_d3_oform'] + OPT[1:]:
    p = RES / c / 'opt.xyz'
    lines = p.read_text().splitlines()
    L += [f'**{STRUCT[c]}** ({lines[0].strip()} atoms; Å; result directory `{c}`)', '', '```']
    L += [l for l in lines[2:2 + int(lines[0])]] + ['```', '']
(ST / 'coordinates_S2.md').write_text('\n'.join(L))

# S3 full bond comparison
MS = [('pcm', 'IEF-PCM'), ('gas', 'Gas'), ('pcm_b3lyp', 'IEF-PCM, no D3'), ('pcm_tz', 'IEF-PCM, def2-TZVP'), ('periodic_i', 'Periodic')]
rows = []
for b in D3['bonds']:
    # X-ray column: the full-precision value of facts.json where the deposited length is a rounding tie (P1-O6,
    # 1.4975(14) in the CIF), so that the table prints what the text prints (review of 5 Oct 2026)
    fx = F.get(f'xray.{b}', {}).get('value')
    xr = fx if isinstance(fx, float) and abs(round(fx, 4) - fx) > 1e-9 else XR[b]
    rows.append([PRETTY(b), f(xr)] + [f(D3['models'][mk]['bond_lengths'][b]) for mk, _ in MS])
table(ST / 'tableS03_bonds.md', ['Bond', 'X-ray (Å)'] + [f'{lab} (Å)' for _, lab in MS],
      [':-----', '---:'] + ['---:'] * len(MS), rows,
      'All bond lengths of the ion pair: X-ray and optimised (B3LYP-D3(BJ)/6-311++G(d,p) unless stated; periodic model '
      'PBE-D3(BJ), formula unit 1).', 's-bonds',
      notes=['In the gas phase H1A is bonded to O5 (1.00 Å), so the N1–H1A entry there is a non-bonded distance. X-ray '
             'C–H riding at 0.93 Å and O–H restrained to 0.82 Å.'], doc='si', section='S3')

# S4 angles and torsions
rows = []
for a in D3['models']['pcm']['angle_devs']:
    xs = None
    rows.append([a.replace('-', '–')] + [f(D3['models'][mk]['angle_devs'][a], 1, sign=True) for mk, _ in MS])
tors = [('C3–C4–S1–N2', 'tau_C3_C4_S1_N2')]
table(ST / 'tableS04_angles.md', ['Angle', ] + [f'{lab} − X-ray (°)' for _, lab in MS],
      [':-----'] + ['---:'] * len(MS), rows,
      'Deviations of the 22 heavy-atom bond angles from the X-ray structure (model − X-ray).', 's-angles',
      notes=[f"C3–C4–S1–N2 torsion: X-ray {f(v('xray.tau_C3C4S1N2'), 1)}°, IEF-PCM {f(v('pcm.tau_C3C4S1N2'), 1)}°, gas "
             f"{f(v('gas.tau_C3C4S1N2'), 1)}°, no D3 {f(v('sens.b3lyp.tau_C3C4S1N2'), 1)}°, def2-TZVP {f(v('sens.tz.tau_C3C4S1N2'), 1)}°."], doc='si', section='S3')

# S5 frequencies
fr = J(RES / 'sfa_pair_pcm_d3/freq.json')
SF = v('freq.scale_factor')
rows = []
for m in D5['modes']:
    ch = ', '.join(f"{a.replace('str ', 'ν ').replace('bend ', 'δ ').replace('-', '–').replace('...', '···')} {w:.2f}" for a, w in m['top_internal'][:3])
    fj = fr['ir'][m['n'] - 1]
    rows.append([m['n'], f(fj['freq_cm-1'], 0), f(fj['freq_cm-1'] * SF, 0), f(fj['intensity_km_mol'], 1), ch])
table(ST / 'tableS05_freq.md', ['Mode', '*ν* (cm⁻¹)', 'Scaled *ν* (cm⁻¹)', 'IR intensity (km mol⁻¹)', 'Largest internal-coordinate contributions'],
      ['---:', '---:', '---:', '---:', ':----'], rows,
      'Harmonic frequencies and IR intensities of the ion pair in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p), finite-'
      'difference Hessian (step 0.005 Bohr), no imaginary frequency. Scaled with 0.9688, the fundamental-frequency factor '
      'of Merrick et al. [@merrick2007evaluation] for B3LYP/6-311+G(d,p), the closest basis in their table (6-311++G(d,p) adds only '
      'diffuse s functions on H). Contributions are squared projections of each mode on bond stretches (ν) and angle bends (δ), '
      'normalised to 1 (a heuristic character, not a potential-energy distribution).', 's-freq',
      notes=['Mode 65 (2190 cm⁻¹) is the N1–H1A stretch of the N1–H1A···O5 bond (97 % N1–H1A stretch character); '
             'the value is harmonic, the mode strongly anharmonic. The O–H stretches (modes 74–75) belong to O–H '  # prose audit b3r2: wording
             'groups that are free in the ion pair but hydrogen-bonded in the crystal.'], doc='si', section='S4')

# S6 thermochemistry
th = D5['thermo']
lo_s, lo_n = v('freq.lowest3'), v('pcm.oform.lowest3')
rows = [['*d*(N1–H1A) / *d*(H1A···O5) (Å)', f"{f(v('proton.pcm.d_N1_H1A'))} / {f(v('proton.pcm.d_H1A_O5'))}",
         f"{f(v('proton.scanO.d_N1_H1A'))} / {f(v('proton.scanO.d_H1A_O5'))}", '–'],
        ['*E*~el~ (Ha)', f(v('thermo.E'), 6), f(v('pcm.oform.E'), 6), f(v('pt.pcm.dE'), 2)],
        ['Imaginary frequencies', v('freq.n_imag'), v('pcm.oform.n_imag'), '–'],
        ['Three lowest frequencies (cm⁻¹)', ', '.join(f(x_, 1) for x_ in lo_s), ', '.join(f(x_, 1) for x_ in lo_n), '–'],
        ['Proton stretch (cm⁻¹; scaled; IR intensity in km mol⁻¹)', f"N1–H1A {f(v('freq.NH_Hbond'), 0)}; {f(v('freq.NH_Hbond.scaled'), 0)}; {f(v('freq.NH_Hbond.I'), 0)}",
         f"O5–H1A {f(v('pcm.oform.nu_OH'), 0)}; {f(v('pcm.oform.nu_OH.scaled'), 0)}; {f(v('pcm.oform.I_OH'), 0)}", f(v('pt.pcm.nu_shift'), 0, sign=True)],
        ['ZPE (Ha)', f(v('thermo.ZPE_Ha'), 6), f(v('pcm.oform.ZPE_Ha'), 6), f(v('pt.pcm.dZPE'), 2)],
        ['ZPE (kcal mol⁻¹)', f(v('thermo.ZPE'), 2), f(v('pcm.oform.ZPE'), 2), f(v('pt.pcm.dZPE'), 2)],
        ['*E*~el~ + ZPE', '–', '–', f(v('pt.pcm.dEZPE'), 2)],
        ['*H*(298.15 K) (Ha)', f(v('thermo.H298'), 6), f(v('pcm.oform.H298'), 6), f(v('pt.pcm.dH298'), 2)],
        ['*G*(298.15 K) (Ha)', f(v('thermo.G298'), 6), f(v('pcm.oform.G298'), 6), f(v('pt.pcm.dG298'), 2)],
        ['Thermal correction to *H* (kcal mol⁻¹)', f(v('thermo.Hcorr'), 2), f(v('pcm.oform.Hcorr'), 2), '–'],
        ['Thermal correction to *G* (kcal mol⁻¹)', f(v('thermo.Gcorr'), 2), f(v('pcm.oform.Gcorr'), 2), '–'],
        ['*S* (J mol⁻¹ K⁻¹)', f(v('thermo.S'), 1), f(v('pcm.oform.S'), 1), f(v('pt.pcm.dS'), 1, sign=True)],
        ['*S*~trans~ / *S*~rot~ / *S*~vib~ (J mol⁻¹ K⁻¹)', f"{f(th['S_tr'], 1)} / {f(th['S_rot'], 1)} / {f(th['S_vib'], 1)}", '–', '–']]
table(ST / 'tableS06_thermo.md', ['Quantity', 'Salt minimum', 'Neutral minimum', 'Neutral − salt^a^'], [':---------', '---:', '---:', '---:'], rows,
      'Harmonic frequencies and thermochemistry of the two minima of the ion pair in IEF-PCM (water), '
      'B3LYP-D3(BJ)/6-311++G(d,p): the salt and the neutral, O5-protonated pair at the end of the relaxed scan. Finite-difference '
      'Hessians (162 gradients each), rigid-rotor, harmonic-oscillator, ideal-gas partition functions at 298.15 K and 1 atm on '
      'top of the PCM electronic energy, unscaled frequencies, no low-mode correction.', 's-thermo',
      notes=['^a^ kcal mol⁻¹, except the proton-stretch row (cm⁻¹, unscaled) and *S* (J mol⁻¹ K⁻¹). Scaled wavenumbers use 0.9688 '
             '[@merrick2007evaluation]. The proton stretch is the strongest IR band of each minimum; the O–H stretches of the anion '
             f"half lie at {f(v('freq.OH_1'), 0)} and {f(v('freq.OH_2'), 0)} cm⁻¹ in the salt and at {' and '.join(f(x_, 0) for x_ in v('pcm.oform.anion_OH'))} cm⁻¹ "
             'in the neutral pair. Both proton stretches are harmonic values of strongly anharmonic modes, but the sign of ΔZPE does not '
             f"rest on them: their zero-point energies favour the salt by {f(v('pt.pcm.dZPE_proton'), 2)} kcal mol⁻¹, and without them ΔZPE is "
             f"{f(v('pt.pcm.dZPE_wo_proton'), 2)} kcal mol⁻¹."],  # prose audit b3r2: re-review O8 (facts pt.pcm.dZPE_proton, _wo_proton)
      doc='si', section='S4')

# S7 scan
rows = []
for p in D4['points']:
    rows.append([p['chain'], p['index'], p['kind'], f(p['target'], 2) if p['target'] else '–', f(p['d_N1_H1A']), f(p['d_H1A_O5']),
                 f(p['delta'], 3, sign=True), f(p['d_N1_O5']), f(p['E_rel_kcal'], 2), 'yes' if p['converged'] else 'no^a^', p['steps']])
table(ST / 'tableS07_scan.md', ['Chain', 'Point', 'Kind', 'Target (Å)', '*d*(N1–H1A) (Å)', '*d*(H1A···O5) (Å)', '*δ* (Å)', '*d*(N1···O5) (Å)',
                                '*E*~rel~ (kcal mol⁻¹)', 'Converged', 'Steps'],
      [':---', '---:', ':---', '---:', '---:', '---:', '---:', '---:', '---:', ':---', '---:'], rows,
      'Relaxed proton scan along *d*(N1–H1A) in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p): chain A from the salt '
      'minimum, chain B from the neutral minimum, 0.10 Å steps. Electronic energies relative to the salt minimum (A, 0).', 's-scan',
      notes=[f"^a^ 61 steps; the constrained-projected gradient at the last step (max {f(v('scan.A02_gmax_proj') * 1e4, 1)} × 10⁻⁴ Ha Bohr⁻¹) is below the "
             f"4.5 × 10⁻⁴ threshold and the last step changed *E* by {f(v('scan.A02_last_dE'), 4)} kcal mol⁻¹. At 1.35 Å, the only point reached "
             f"from both ends, the chains differ by {f(v('scan.hysteresis_1.35'), 3)} kcal mol⁻¹. Barrier from the salt side "
             f"{f(v('scan.barrier'), 2)} kcal mol⁻¹ (cubic spline {f(v('scan.barrier_spline'), 2)}), from the neutral side {f(v('scan.reverse_barrier'), 2)} kcal mol⁻¹."],
      doc='si', section='S5')

# S8 dielectric ladder
rows = []
for eps, k, src in [('1 (gas)', 'gas', 'X-ray'), ('4', 'eps4', 'water optimum'), ('8', 'eps8', 'water optimum'), ('20', 'eps20', 'water optimum'), ('78.4 (water)', 'pcm', 'X-ray')]:
    rows.append([eps, src, f(v(f'proton.{k}.d_N1_H1A')), f(v(f'proton.{k}.d_H1A_O5')), f(v(f'proton.{k}.delta'), 3, sign=True),
                 f(v(f'proton.{k}.d_N1_O5')), f(v(f'proton.{k}.P1O5')), f(v(f'proton.{k}.C1N1')), 'N1' if v(f'proton.{k}.side') == 'N' else 'O5'])
table(ST / 'tableS08_eps.md', ['*ε*', 'Start', '*d*(N1–H1A) (Å)', '*d*(H1A···O5) (Å)', '*δ* (Å)', '*d*(N1···O5) (Å)', '*d*(P1–O5) (Å)', '*d*(C1–N1) (Å)', 'Proton on'],
      [':---', ':---', '---:', '---:', '---:', '---:', '---:', '---:', ':---'], rows,
      'Dielectric ladder: the ion pair optimised in IEF-PCM with *ε* = 1–78.4, B3LYP-D3(BJ)/6-311++G(d,p). Heavy-atom '
      'markers independent of the H position: P1–O5 lengthens and C1–N1 shortens when O5 is protonated.', 's-eps',
      notes=['The salt minimum disappears between *ε* = 4 and 8. Energies at different *ε* are not comparable and are not listed.'], doc='si', section='S5')

# S9 functional check (from log)
rows = []
for lab, mk in [('B3LYP-D3(BJ)', 'b3lyp'), ('ωB97X-D', 'wb97xd'), ('PBE0-D3(BJ)', 'pbe0'), ('M06-2X', 'm062x'), ('DF-MP2', 'mp2')]:
    rows.append([lab, f(v(f'fcheck.{mk}.top_N'), 2, sign=True), f(v(f'fcheck.{mk}.dE_O_N'), 2, sign=True), f(v(f'fcheck.{mk}.top_O'), 2, sign=True)])
table(ST / 'tableS09_fcheck.md', ['Method', '*E*(1.25 Å) − *E*(salt)', '*E*(neutral) − *E*(salt)', '*E*(1.25 Å) − *E*(neutral)'],
      [':------', '---:', '---:', '---:'], rows,
      'Functional check (kcal mol⁻¹): single points in IEF-PCM (water) with 6-311++G(d,p) on the B3LYP-D3(BJ) geometries of '
      'the salt minimum (*d*(N1–H1A) = 1.098 Å), the scan point at 1.25 Å and the neutral minimum (1.631 Å).', 's-fcheck',
      notes=['Electronic energies, not re-optimised with each method. ωB97X-D includes its Chai–Head-Gordon dispersion term; DF-MP2 '
             # refs audit round 2 (29 Sep 2026): cite the DF-MP2 method and its fitting basis here too (refs_round1 minor 10)
             '[@feyereisen1993use] is frozen-core (aug-cc-pVTZ-RI fitting basis [@weigend2002efficient]) with the IEF-PCM reaction '
             'field of the RHF density. All five methods put the '
             'neutral pair below the salt. A negative first column means that on these geometries the salt is not separated from '
             'the neutral pair by a barrier; the largest value, ' + f(v('fcheck.wb97xd.top_N'), 2) + ' kcal mol⁻¹ (ωB97X-D), is below the '
             + f(v('zpe.NH'), 2) + ' kcal mol⁻¹ zero-point energy of the N1–H1A stretch. The four DFT sets reproduce an earlier run on '
             'the same geometries within ' + f(v('fcheck.vs_burst1_max_Ha') * 1e11, 1) + ' × 10⁻¹¹ Ha.'], doc='si', section='S5')

# S10 ZPE vs barrier
rows = [['Electronic barrier, salt → neutral (scan grid)', f(v('scan.barrier'), 2)], ['Electronic barrier, cubic spline', f(v('scan.barrier_spline'), 2)],
        ['*E*(neutral) − *E*(salt)', f(v('scan.Omin'), 2)], ['½*hν*, N1–H1A stretch normal mode (2190 cm⁻¹)', f(v('zpe.NH'), 2)],
        ['½*hν*, N1–H1A adiabatic local mode (1896 cm⁻¹)', f(v('zpe.local_mode_zpe'), 2)],
        ['Frequency from the scan curvature (cm⁻¹)', f(v('zpe.scan_curvature'), 0)],
        ['1-D ground level above the salt minimum', f(v('zpe.oned_E0'), 2)], ['1-D ground level minus barrier top', f(v('zpe.oned_E0_minus_top'), 2)],
        ['1-D ground-state probability on the O5 side', f(v('zpe.oned_P_O'), 2)]]
table(ST / 'tableS10_zpe.md', ['Quantity', 'Value (kcal mol⁻¹ unless stated)'], [':---------', '---:'], rows,
      'Zero-point energy of the N1–H1A stretch against the proton-transfer barrier in IEF-PCM (water), B3LYP-D3(BJ)/'
      '6-311++G(d,p). The 1-D model solves the proton motion on the relaxed scan profile (effective mass 0.940 u, '
      'sinc-DVR).', 's-zpe',
      notes=['Harmonic, one-dimensional estimates; the other modes are not followed along the path.'], doc='si', section='S5')

# S11 intra-pair contacts in every model + sensitivity
rows = []
hbp = {h['name']: h for h in D8['structures']['relax_i']['hbonds']}
for lab, k, mk in [('X-ray', 'xray', 'xray'), ('Gas phase, B3LYP-D3(BJ)/6-311++G(d,p)', 'gas', 'gas'),
                   ('IEF-PCM, B3LYP-D3(BJ)/6-311++G(d,p)', 'pcm', 'pcm'), ('IEF-PCM, B3LYP/6-311++G(d,p) (no D3)', 'sens.b3lyp', 'pcm_b3lyp'),
                   ('IEF-PCM, B3LYP-D3(BJ)/def2-TZVP', 'sens.tz', 'pcm_tz'), ('Periodic model, PBE-D3(BJ)', 'periodic', 'periodic_i')]:
    g = D3['models'][mk]
    pr = g['proton']
    if mk == 'periodic_i':
        a1, a2 = hbp['N1-H1A...O5']['DHA'], hbp['C6-H6...O4']['DHA']
        tors = f(pr['tau_C3_C4_S1_N2'], 1)
    else:
        a1, a2 = pr['ang_N1_H1A_O5'], pr['ang_C6_H6_O4']
        tors = f(pr['tau_C3_C4_S1_N2'], 1)
    stats_ = ['–', '–', '–'] if mk == 'xray' else [f(g['bond_stats']['_heavy']['mae']), f(g['bond_stats']['_heavy']['rmsd']), f(g['angle_stats_heavy']['mae'], 1)]
    # X-ray δ to 0.01 Å, as in Table 2 (numbers audit round 1): its third decimal is not supported by the H position.
    rows.append([lab, f(v(f'proton.{k}.delta'), 2 if k == 'xray' else 3, sign=True), f(v(f'proton.{k}.d_N1_O5')), f(a1, 1), f(pr['d_C6_O4']), f(a2, 1)] + stats_ + [tors])
table(ST / 'tableS11_sensitivity.md', ['Model', '*δ* (Å)', '*d*(N1···O5) (Å)', '∠N1–H1A···O5 (°)', '*d*(C6···O4) (Å)', '∠C6–H6···O4 (°)',
                                       'Heavy-bond MAE (Å)', 'Heavy-bond RMSD (Å)', 'Angle MAE (°)', 'C3–C4–S1–N2 (°)'],
      [':--------------------------', '--------------:', '------------------:', '-------------------------:', '------------------:', '---------------------:', '--------------:', '--------------:', '-------------:', '-------------------------:'], rows,
      'The two contacts that close the R²₂(8) ring of the ion pair, and the agreement of every model with the X-ray '
      'geometry. The models differ in environment (gas phase, IEF-PCM, crystal), dispersion (with and without D3) and '
      'basis set (6-311++G(d,p), def2-TZVP).', 's-sensitivity',
      notes=['In the gas phase H1A sits on O5, so the first contact is O5–H1A···N1 there. Without D3 the C6–H6···O4 '
             'contact opens and the ring breaks while the proton position is unchanged. Only the periodic model keeps the '
             'crystal torsion; every isolated ion pair relaxes it to about 90°.'], doc='si', section='S6')

# S12 periodic convergence
conv = D8['convergence']
rows = [[f"{c_['ecut_Ry'][0]}/{c_['ecut_Ry'][1]}", c_['k'].replace('x', '×'), f(v(f'bulk.qeconv.{i}.dE'), 3), f(c_['max_dF_Ry_bohr'], 6), f(c_['dP_kbar'], 2)]
        for i, c_ in enumerate(conv['ecut'] + conv['kgrid'])]
table(ST / 'tableS12_periodic_conv.md', ['Cutoff (Ry)', '*k* grid', 'Δ*E* (meV atom⁻¹)', 'max Δ*F* (Ry Bohr⁻¹)', 'Δ*P* (kbar)'],
      [':---', ':---', '---:', '---:', '---:'], rows,
      'Convergence of the periodic PBE-D3(BJ) model at the X-ray geometry: cutoff at a 3×2×2 grid (first three rows, '
      'reference 70/560 Ry) and *k* grid at 50/400 Ry (last four rows, reference 5×2×3). Chosen: 70/560 Ry and 2×1×1.',
      's-periodic-conv',
      notes=['The chosen cutoff is the highest one tested, so the table does not show the total energy converged in the '
             'cutoff; what the paper uses is an energy difference, and a second pseudopotential set at another cutoff '
             'changes it by at most 0.8 kJ mol⁻¹ per formula unit (Section 4.3 of the main text).'], doc='si', section='S7')

# S13 periodic relaxations
pq = J(RES / 'periodic_qe/periodic_qe.json')
rows = []
PBE0_KJ = {'i': '0 (reference)', 'ii': f(v('periodic.hybrid.ii_best_kJ'), 2, sign=True), 'iirfix': f(v('periodic.hybrid.e100.best_kJ'), 2, sign=True)}
for k, lab in [('i', 'From the X-ray structure (i)'), ('ii', 'From H1A on O5 in both formula units (ii)'),
               ('iirfix', 'Constrained single configuration: H1A and O5 held^a^')]:  # prose audit burst 3: run name dropped
    s_ = D8['structures'][f'relax_{k}']
    rows.append([lab, f(pq['relax'][k]['energy_Ry'], 6), ('0 (reference)' if k == 'i' else f((pq['relax'][k]['energy_Ry'] - pq['relax']['i']['energy_Ry']) * 13.605693 * 96.485332 / 2, 2, sign=True)),
                 PBE0_KJ[k], f(s_['fu1']['delta'], 3, sign=True), f(s_['fu1']['d_N1_O5']), f(s_['rmsd_heavy_shift_removed']), pq['relax'][k]['bfgs_steps'][1],
                 f(v('periodic.iirfix.max_force_free'), 5) + '^b^' if k == 'iirfix' else f(pq['relax'][k]['max_force_Ry_bohr'], 5)])
    if k == 'ii':
        # single-transfer control (VM burst 3): H1A on O5 in formula unit 1 only, unconstrained BFGS without symmetry
        # prose audit burst 3: label (iii) as in main 2.3, not the run name iir1
        rows.append(['From H1A on O5 in formula unit 1 only (iii)^c^', f(v('periodic.iir1.E'), 6), f(v('periodic.iir1.dE_kJ'), 2, sign=True),
                     '–', f(v('periodic.iir1.delta'), 3, sign=True), f(v('periodic.iir1.d_N1_O5')), f(v('periodic.iir1.rmsd_heavy')), v('periodic.iir1.steps'),
                     f(v('periodic.iir1.max_force'), 5)])
table(ST / 'tableS13_periodic_relax.md', ['Start', '*E* (Ry)', 'Δ*E*, PBE-D3(BJ) (kJ mol⁻¹ per formula unit)', 'Δ*E*, PBE0-D3(BJ) (kJ mol⁻¹ per formula unit)',
                                          '*δ* (Å)', '*d*(N1···O5) (Å)', 'Heavy-atom RMSD to X-ray (Å)', 'BFGS steps', 'Max force (Ry Bohr⁻¹)'],
      [':-------------------------------', '----------------------------:', '-----------------------:', '-----------------------:', '---------------:', '-----------------------:', '---------------:', '-------------:', '--------------------:'], rows,
      'Periodic relaxations of the *Pc* cell and their energies relative to (i). PBE-D3(BJ) with the cell fixed at the '
      'X-ray values, 70/560 Ry, 2×1×1. The PBE0-D3(BJ) correction is the difference of PBE0-D3(BJ) and PBE-D3(BJ) single '
      'points at the relaxed geometries (SG15 potentials, 100/400 Ry), added to the PBE-D3(BJ) difference.', 's-periodic-relax',
      notes=[f"^a^ H1A of both formula units placed {f(v('protocol.qe_start_OH'), 2)} Å from O5 in the relaxed cell (i), then H1A and O5 "
             f"held in space while every other atom relaxed: a constrained structure, not a minimum; the lattice moved by "
             f"{f(v('periodic.iirfix.shift'), 2)} Å only. ^b^ Largest force on the free atoms; the held H1A carries "
             f"{f(v('periodic.iirfix.F_H1A'), 3)} Ry Bohr⁻¹ ({f(v('periodic.iirfix.F_H1A_toN1_eVA'), 1)} eV Å⁻¹ towards N1). The hybrid correction is "
             f"{f(v('periodic.hybrid.e100.corr_kJ'), 1)} kJ mol⁻¹ per formula unit for the constrained neutral "
             f"({f(v('periodic.hybrid.e80.corr_kJ'), 1)} at 80/320 Ry, {f(v('periodic.hybrid.nq211.corr_kJ'), 1)} with a 2×1×1 *q* mesh), "
             f"{f(v('periodic.hybrid.e100.corr_exx_meV'), 1)} of its {f(v('periodic.hybrid.e100.corr_meV'), 1)} meV per cell from exact "
             f"exchange. Start (ii) relaxes to the same structure as (i). A first constrained run that held H1A alone let the lattice "
             f"slide by {f(v('periodic.iifix_shift_mag'), 2)} Å until N1 bonded to H1A again; it is not used. The single configuration "
             f"lies {f(v('periodic.profile.iirfix_excess'), 1)} kJ mol⁻¹ above the constrained PBE-D3(BJ) N1–H1A profile at its N1–H1A distance "
             f"(@tbl:s-periodic-profile), which replaces it as the measure of the cost. ^c^ No symmetry imposed ({v('periodic.iir1.n_symops')} symmetry operation); "
             'the proton of formula unit 1 returns to N1 and formula unit 2 stays a salt, so the cell relaxes to the salt of (i) '
             f"(all-atom RMSD {f(v('periodic.iir1.rmsd_vs_i'), 3)} Å from it). No hybrid single point was run for this cell."], doc='si', section='S7')

# S13b constrained periodic N1-H1A profile (referee M2/M3; facts section 32) -- printed as Table S15
PK = ['p110', 'p125', 'p140', 'p155', 'p170']
HBN = [('N1H1B_O6', 'N1–H1B···O6'), ('N1H1C_O6', 'N1–H1C···O6'), ('N2H2A_O3', 'N2–H2A···O3'), ('N2H2B_O2', 'N2–H2B···O2'),
       ('O3H3A_O1', 'O3–H3A···O1'), ('O4H4A_O5', 'O4–H4A···O5'), ('C3H3_O6', 'C3–H3···O6'), ('C6H6_O4', 'C6–H6···O4')]
rows = [['*d*(N1–H1A) (Å)', f(v('proton.periodic.d_N1_H1A'))] + [f(v(f'proton.periodic.{k}.d_N1_H1A')) for k in PK],
        ['*d*(H1A···O5) (Å)', f(v('proton.periodic.d_H1A_O5'))] + [f(v(f'proton.periodic.{k}.d_H1A_O5')) for k in PK],
        ['*δ* (Å)', f(v('proton.periodic.delta'), 3, sign=True)] + [f(v(f'proton.periodic.{k}.delta'), 3, sign=True) for k in PK],
        ['*d*(N1···O5) (Å)', f(v('proton.periodic.d_N1_O5'))] + [f(v(f'proton.periodic.{k}.d_N1_O5')) for k in PK],
        ['*d*(P1–O5) (Å)', f(v('proton.periodic.P1O5'))] + [f(v(f'proton.periodic.{k}.P1O5')) for k in PK]]
for hk, hl in HBN:
    rows.append([f'*d*(D···A) {hl} (Å)', f(v(f'periodic.hb.{hk}.DA'))] + [f(v(f'periodic.profile.{k}.hb.{hk}.DA')) for k in PK])
rows += [['Δ*E*, PBE-D3(BJ) (kJ mol⁻¹)', '0'] + [f(v(f'periodic.profile.{k}.PBE_kJ'), 1, sign=True) for k in PK],
         ['Hybrid correction (kJ mol⁻¹)', '–'] + [f(v(f'periodic.profile.{k}.corr_kJ'), 1, sign=True) for k in PK],
         ['Δ*E*, PBE0-D3(BJ) (kJ mol⁻¹)', '0'] + [f(v(f'periodic.profile.{k}.PBE0_kJ'), 1, sign=True) for k in PK],
         ['d*E*/d*d*(N1–H1A) (kJ mol⁻¹ Å⁻¹)', '–'] + [f(v(f'periodic.profile.{k}.slope_kJ'), 0, sign=True) for k in PK],
         ['Residual pressure (kbar)', f(v('periodic.salt.pressure'), 1)] + [f(v(f'periodic.profile.{k}.pressure'), 1, sign=True) for k in PK],
         ['Heavy-atom RMSD to X-ray (Å)', f(v('periodic.rmsd_heavy'))] + [f(v(f'periodic.profile.{k}.rmsd_heavy')) for k in PK],
         ['Optimisation steps', v('periodic.bfgs_steps_i')] + [v(f'periodic.profile.{k}.steps') for k in PK]]
table(ST / 'tableS13b_periodic_profile.md', ['Quantity', 'Salt (i)'] + [f"{f(t_, 2)} Å" for t_ in v('periodic.profile.targets')],
      [':---', '---:', '---:', '---:', '---:', '---:', '---:'], rows,
      'Constrained proton-transfer profile in the crystal: *d*(N1–H1A) of both formula units held at 1.10–1.70 Å by a '
      'distance constraint, every other atom relaxed (PBE-D3(BJ), SSSP, 70/560 Ry, 2×1×1, cell fixed at the X-ray values, '
      'damped dynamics), with PBE0-D3(BJ) single points at the same geometries. Energies per formula unit relative to the '
      'salt (i); geometry of formula unit 1 (formula unit 2 is identical).', 's-periodic-profile',
      notes=['Hybrid correction: *E*(PBE0-D3(BJ)) − *E*(PBE-D3(BJ)) differences from SG15 single points (100/400 Ry, *q* 1×1×1), added to '
             f"the SSSP PBE-D3(BJ) difference; the D3(BJ) parameters contribute at most {f(v('periodic.profile.corr_D3_max'), 1)} meV per cell, the rest is "
             f"exact exchange. d*E*/d*d* is the constraint force at the constrained minimum, (*F*~N1~ − *F*~H1A~)·*u*/2 per formula unit: positive "
             'at every point, so the held proton is pushed back towards N1 and no point is a minimum. Converged on the constrained forces '
             f"below {f(v('periodic.profile.forc_conv_thr'), 3)} Ry Bohr⁻¹. At the "
             f"N1–H1A distance of the single configuration with O5–H1A held at 0.98 Å ({f(v('proton.periodic.iirfix.d_N1_H1A'), 3)} Å, @tbl:s-periodic-relax) the PBE-D3(BJ) profile "
             f"interpolates to {f(v('periodic.profile.at_iirfix_d'), 1, sign=True)} kJ mol⁻¹, {f(v('periodic.profile.iirfix_excess'), 1)} kJ mol⁻¹ below that "
             'configuration. The cell was not relaxed; judged from the residual pressures with a bulk modulus of 10 GPa (an estimate, not a '
             f"calculation), relaxing it would lower the salt by about {f(abs(v('periodic.profile.stress_B10')), 1)} kJ mol⁻¹ per formula unit more than "
             'the 1.70 Å cell and so raise the cost slightly.'],
      doc='si', section='S7')

# S14 periodic H-bonds
rows = []
for h in D8['structures']['relax_i']['hbonds']:
    x = hb_x[h['name']]
    rows.append([h['name'].replace('...', '···').replace('-', '–'), SYMTXT[h['sym']], x['cif']['DA'], f(h['DA']), f(h['dDA'], 3, sign=True),
                 x['cif']['HA'], f(h['HA']), x['cif']['DHA'], f(h['DHA'], 1)])
table(ST / 'tableS14_periodic_hbonds.md', ['Contact', 'Acceptor at', '*d*(D···A) X-ray (Å)', '*d*(D···A) periodic (Å)', 'Δ (Å)', '*d*(H···A) X-ray (Å)', '*d*(H···A) periodic (Å)', 'Angle X-ray (°)', 'Angle periodic (°)'],
      [':--------------------------------', ':------------------------------------', '-------------------:', '-------------------:', '----------------:', '-------------------:', '-------------------:', '---------------:', '---------------:'], rows,
      'The nine hydrogen bonds of the crystal: deposited values and the relaxed periodic PBE-D3(BJ) model (i).', 's-periodic-hbonds',
      notes=[f"O4–H4A···O5: the periodic model puts H4A at {f(v('periodic.O4H4A'))} Å from O4 and {f(v('periodic.H4AO5'))} Å from O5 "
             f"(fraction {f(v('periodic.O4H4A_frac'), 2)} of the O4···O5 line), against the restrained X-ray 0.80/1.76 Å."], doc='si', section='S7')

# S15 descriptors, DSCF vs Koopmans
rows = []
for sp, mod, dk in [('Ion pair', 'IEF-PCM (water)', 'pcm'), ('Cation', 'IEF-PCM (water)', 'cation.pcm'), ('Cation', 'Gas phase', 'cation.gas'),
                    ('Sulfanilamide', 'IEF-PCM (water)', 'sfa.pcm'), ('Sulfanilamide', 'Gas phase', 'sfa.gas')]:
    for meth, pre in [('ΔSCF', f'dscf.{dk}.'), ('Koopmans', f'dscf.{dk}.koop.')]:
        rows.append([sp if meth == 'ΔSCF' else '', mod if meth == 'ΔSCF' else '', meth, f(v(pre + 'I')), f(v(pre + 'A')), f(v(pre + 'eta')),
                     f(v(pre + 'chi')), f(v(pre + 'omega')), f(v(pre + 'S'))])
table(ST / 'tableS15_descriptors.md', ['Species', 'Model', 'Method', '*I* (eV)', '*A* (eV)', '*η* (eV)', '*χ* (eV)', '*ω* (eV)', '*S* (eV⁻¹)'],
      [':---', ':---', ':---', '---:', '---:', '---:', '---:', '---:', '---:'], rows,
      'Global reactivity descriptors from vertical ΔSCF energies and from the orbital energies (Koopmans-type, *I* = −*E*~HOMO~, '
      '*A* = −*E*~LUMO~), B3LYP-D3(BJ)/6-311++G(d,p).', 's-descriptors',
      notes=['Condensed Fukui functions of sulfanilamide and the cation are in @tbl:s-fukui.'], doc='si', section='S9')

# S16 atomic charges
ch = J(RES / 'sfa_pair_pcm_d3/mwfn_extra/charges.json')
labs = [l.split(',')[1].replace('#1_554', '') for l in (RES / 'sfa_pair_pcm_d3/mwfn_extra/charges.csv').read_text().splitlines()[1:]]
rows = [[a] + [f(ch['atomic'][s_][i]) for s_ in ['NPA', 'Hirshfeld', 'ADCH', 'CM5']] for i, a in enumerate(labs)]
fs = ch['fragment_sums']
rows.append(['Σ cation'] + [f(fs['4-sulfamoylanilinium cation'][s_]) for s_ in ['NPA', 'Hirshfeld', 'ADCH', 'CM5']])
table(ST / 'tableS16_charges.md', ['Atom', 'NPA (e)', 'Hirshfeld (e)', 'ADCH (e)', 'CM5 (e)'], [':---', '---:', '---:', '---:', '---:'], rows,
      'Atomic charges of the ion pair in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p). NPA from PySCF natural atomic '
      'orbitals (not NBO 7); Hirshfeld, ADCH and CM5 from Multiwfn.', 's-charges',
      notes=['The ADCH values of C3, C6, H3 and H6 have the wrong sign for a C–H bond and should not be interpreted '
             'atom by atom.'], doc='si', section='S9')

# S17 MEP extrema
ms = J(RES / 'sfa_pair_pcm_d3/mwfn_extra/mep_summary.json')
rows = [[e['type'] + (' (global)' if e['global'] == 'yes' else ''), f(e['ESP_kcal_mol'], 1), f(e['ESP_au'], 4), e['nearest_atoms'].replace('#1_554', '').replace(' A)', ' Å)')]
        for e in sorted(ms['extrema'], key=lambda e: e['ESP_kcal_mol'])]
table(ST / 'tableS17_mep.md', ['Extremum', '*V* (kcal mol⁻¹)', '*V* (a.u.)', 'Nearest atoms (Å)'], [':---', '---:', '---:', ':---'], rows,
      'Extrema of the electrostatic potential on the 0.001 a.u. electron-density surface of the ion pair in IEF-PCM '
      f"(water), B3LYP-D3(BJ)/6-311++G(d,p). Surface area {f(v('mep.pcm.area_A2'), 1)} Å², polar fraction {f(v('mep.pcm.polar_area_pct'), 1)} %, "
      f"*σ*²~tot~ = {f(v('mep.pcm.sigma2_tot_kcal2'), 0)} (kcal mol⁻¹)², *ν* = {f(v('mep.pcm.nu'), 3)}, molecular polarity index {f(v('mep.pcm.MPI_eV'), 2)} eV.",
      's-mep', notes=['Two of the listed "maxima" are negative-valued saddle-like points on the anion side, as reported by Multiwfn.'], doc='si', section='S9')

# S18 QTAIM and IGMH
rows = []
for c, lab in [('sfa_pair_pcm_d3', 'IEF-PCM'), ('sfa_pair_gas_d3', 'Gas phase')]:
    for h in D6['qtaim'][c]['hbonds']:
        nm = h['name'].split(' (')[0].replace('#1_554', '').replace('...', '···').replace('-', '–')
        rows.append([lab, nm, f(h['DA']), f(h['DHA'], 1), f(h['rho'], 4), f(h['lap'], 4), f(h['Hd'], 4), f(h['mGV'], 2),
                     f(h['E_HB_kcal'], 1), f(h['BE_em_neutral'], 1), f(h['BE_em_charged'], 1)])
table(ST / 'tableS18_qtaim.md', ['Model', 'Contact', '*d*(D···A) (Å)', 'Angle (°)', '*ρ* (a.u.)', '∇²*ρ* (a.u.)', '*H* (a.u.)', '−*G*/*V*',
                                 '*V*/2 (kcal mol⁻¹)', 'BE neutral fit (kcal mol⁻¹)', 'BE charged fit (kcal mol⁻¹)'],
      [':-----------', ':------------------------', '-----------------:', '-----------:', '------------:', '--------------:', '--------------:', '----------:', '--------------------:', '----------------------------:', '----------------------------:'], rows,
      'QTAIM descriptors at the intermolecular bond critical points of the ion pair, B3LYP-D3(BJ)/6-311++G(d,p), with '
      'three hydrogen-bond energy estimates (Espinosa *V*/2; Emamian fits for neutral and charged complexes).', 's-qtaim',
      notes=[f"Not bond energies: *V*/2 of N1–H1A···O5 in IEF-PCM ({f(v('qtaim.pcm.N1H1A_O5.EHB'), 1)} kcal mol⁻¹) exceeds in magnitude the "
             f"counterpoise-corrected interaction energy of the whole ion pair ({f(v('int.pcm.cp'), 2)} kcal mol⁻¹, @tbl:s-interaction), "
             'dispersion included, so these estimates are not used in the main text. '
             'An ion pair also lies outside both Emamian training sets. In the gas phase the proton is on O5, '
             f"so the first contact is O5–H1A···N1. IGMH: integrated *δg*~inter~ {f(v('igmh.pcm.int_dg_inter'), 3)} a.u. (IEF-PCM), "
             f"of which H1A···O5 {f(v('igmh.pcm.H1A_O5'), 1)} %, N1···O5 {f(v('igmh.pcm.N1_O5'), 1)} %, H1A···P1 {f(v('igmh.pcm.H1A_P1'), 1)} %, "
             f"H6···O4 {f(v('igmh.pcm.H6_O4'), 1)} %."], doc='si', section='S9')

# S18b QTAIM of the nine crystal hydrogen bonds in the five-pair cluster (referee M6; facts section 34) -- printed as Table S22
rows = []
for key, h in hb_x.items():
    kq = key.replace('-', '').replace('...', '_')
    nm = key.replace('...', '···').replace('-', '–')
    rng = (f"{f(v(f'qtaim.xclu.{kq}.rho_range_salt')[0], 4)}–{f(v(f'qtaim.xclu.{kq}.rho_range_salt')[1], 4)}"
           if f'qtaim.xclu.{kq}.rho_range_salt' in F else '–')
    rows.append([nm + ('^a^' if kq == 'N1H1A_O5' else ''), f(v(f'qtaim.xclu.{kq}.salt.HA')), f(v(f'qtaim.xclu.{kq}.neutral.HA')),
                 f(v(f'qtaim.xclu.{kq}.salt.rho'), 4), f(v(f'qtaim.xclu.{kq}.neutral.rho'), 4), f(v(f'qtaim.xclu.{kq}.d_rho_pct'), 0, sign=True),
                 f(v(f'qtaim.xclu.{kq}.salt.lap'), 4), f(v(f'qtaim.xclu.{kq}.neutral.lap'), 4),
                 f(v(f'qtaim.xclu.{kq}.salt.H'), 5), f(v(f'qtaim.xclu.{kq}.neutral.H'), 5), v(f'qtaim.xclu.{kq}.n_instances'), rng])
table(ST / 'tableS18b_cluster_qtaim.md', ['Contact', '*d*(H···A), salt (Å)', '*d*(H···A), H1A on O5 (Å)', '*ρ*, salt (a.u.)', '*ρ*, H1A on O5 (a.u.)',
                                          'Δ*ρ* (%)', '∇²*ρ*, salt (a.u.)', '∇²*ρ*, H1A on O5 (a.u.)', '*H*, salt (a.u.)', '*H*, H1A on O5 (a.u.)',
                                          'Instances', '*ρ* over the instances, salt (a.u.)'],
      [':--------------------------', '----------------:', '-----------------:', '------------:', '-------------:', '-------:', '-------------:', '--------------:', '-----------------:', '-----------------:', '----------:', '--------------------------:'], rows,
      'QTAIM descriptors of the nine hydrogen bonds of the crystal in the five-pair cluster. The cluster (135 atoms) was cut '
      'from the relaxed periodic salt and from the periodic cell with H1A held on O5 (O5–H1A 0.98 Å, every pair); '
      'B3LYP-D3(BJ)/6-31+G(d,p) single points in the gas phase at the periodic PBE-D3(BJ) geometries. Values for the contact '
      'at the central ion pair; N1–H1B···O6 has no instance there and is given for one of its two instances, both between outer '
      'pairs.', 's-cluster-qtaim',
      notes=[f"^a^ With H1A on O5 the critical point is that of N1···H1A ({f(v('qtaim.xclu.N1H1A_O5.neutral.HA'), 3)} Å is then O5–H1A). "
             'Instances: occurrences of the contact in the cluster; the last column shows how much the edge of the cluster moves *ρ*. '
             f"The cluster energies are not a model of the crystal: the salt cluster has a HOMO–LUMO gap of {f(v('qtaim.xclu.gap_salt'), 2)} eV "
             f"({f(v('qtaim.xclu.gap_neutral'), 2)} eV with H1A on O5), and in the gas phase the cluster with H1A on O5 lies "
             f"{f(abs(v('qtaim.xclu.dE_kJ')), 1)} kJ mol⁻¹ per ion pair below the salt cluster, whereas the periodic cell puts the same configuration "
             f"{f(v('periodic.iirfix.dE_kJ'), 1, sign=True)} kJ mol⁻¹ per formula unit above the salt."], doc='si', section='S9')

# S19 interaction energies
rows = []
for lab, k in [('IEF-PCM, ion pair (cation + H₂PO₄⁻)', 'int.pcm'), ('Gas phase at the IEF-PCM geometry (cation + H₂PO₄⁻)', 'int.gas_ionic'),
               ('Gas phase, own optimum (sulfanilamide + H₃PO₄)', 'int.gas_neutral')]:
    defs = [kk for kk in F if kk.startswith(k + '.def.')]
    rows.append([lab, f(v(k + '.raw'), 2), f(v(k + '.bsse'), 2), f(v(k + '.cp'), 2), f(v(k + '.d3'), 2),
                 ' + '.join(f(v(d_), 2) for d_ in defs), f(v(k + '.bind'), 2)])
table(ST / 'tableS19_interaction.md', ['Partition', '*E*~int~ (raw)', 'BSSE', '*E*~int~ (CP)', 'of which D3(BJ)', 'Σ*E*~def~', '*E*~bind~'],
      [':---', '---:', '---:', '---:', '---:', '---:', '---:'], rows,
      'Interaction energies of the ion pair (kcal mol⁻¹), B3LYP-D3(BJ)/6-311++G(d,p), counterpoise-corrected; deformation '
      'energies relative to the relaxed fragments in the same environment (cation, then anion; sulfanilamide, then H₃PO₄).', 's-interaction',
      notes=['In IEF-PCM the counterpoise correction is approximate (each species in its own cavity); the gas-phase BSSE at '
             f"the same geometry is {f(v('int.pcm.bsse_gas_same_geom'), 2)} kcal mol⁻¹."], doc='si', section='S9')

# S20 Hirshfeld fingerprint decomposition
hs = J(RES / 'hirshfeld_xray/hirshfeld_summary.json')
cmp_ = J(RES / 'hirshfeld_xray/hirshfeld_compare.json')
keys = []
for ion in ['cation', 'anion']:
    for c_ in hs['targets'][ion]['contacts_reciprocal']:
        if c_['contact'] not in keys:
            keys.append(c_['contact'])
rows = []
for kk in keys:
    cat = {c_['contact']: c_['pct'] for c_ in hs['targets']['cation']['contacts_reciprocal']}
    an = {c_['contact']: c_['pct'] for c_ in hs['targets']['anion']['contacts_reciprocal']}
    rows.append([kk.replace('...', '···'), f(cat.get(kk), 1) if kk in cat else '–', f(an.get(kk), 1) if kk in an else '–'])
HSNOTE = (f"Surface area {f(hs['targets']['cation']['area_A2'], 1)} Å² (cation) and {f(hs['targets']['anion']['area_A2'], 1)} Å² (anion); "
          f"volume {f(hs['targets']['cation']['volume_A3'], 1)} Å³ and {f(hs['targets']['anion']['volume_A3'], 1)} Å³; smallest *d*~norm~ "
          f"{f(hs['targets']['cation']['d_norm_min'], 3)} and {f(hs['targets']['anion']['d_norm_min'], 3)}.")
table(ST / 'tableS20_hirshfeld.md', ['Contact', 'Cation (%)', 'Anion (%)'], [':---', '---:', '---:'], rows,
      'Decomposition of the Hirshfeld surfaces of the cation and the anion into element-pair contacts. Reciprocal contacts '
      'are combined; X-ray structure with neutron-normalised X–H (C–H 1.083, N–H 1.009, O–H 0.983 Å).', 's-hirshfeld',
      notes=[HSNOTE, 'Independent Multiwfn surfaces give the same areas within 1.2 % and every percentage within 0.2 points.'], doc='si', section='S10')

# S21 Hirshfeld red spots
rows = []
for ion in ['cation', 'anion']:
    for rs in hs['targets'][ion]['red_spots']:
        rows.append([ion, f"{rs['inside'].split('#')[0]}···{rs['outside'].split('#')[0]}", f(rs['dist_A']), f(rs['dnorm'], 3)])
table(ST / 'tableS21_hirshfeld_contacts.md', ['Surface', 'Contact (inside···outside)', 'Distance (Å)', '*d*~norm~'], [':---', ':---', '---:', '---:'], rows,
      'Local minima of *d*~norm~ (red spots) on the Hirshfeld surfaces, with the interatomic distance of the contact '
      '(neutron-normalised H positions).', 's-hirshfeld-contacts', doc='si', section='S10')

# S22 TD-DFT ion pair, S23 fragments
for fn, lab, label, cases in [('tableS22_tddft_pair.md', 'ion pair', 's-tddft-pair', [('sfa_pair_pcm_d3/tddft_camb3lyp_noneq', 'CAM-B3LYP'), ('sfa_pair_pcm_d3/tddft_b3lyp_noneq', 'B3LYP')]),
                              ('tableS23_tddft_fragments.md', 'cation and sulfanilamide', 's-tddft-fragments',
                               [('cation_pcm_d3/tddft_camb3lyp_noneq', 'Cation, CAM-B3LYP'), ('sulfanilamide_pcm_d3/tddft_camb3lyp_noneq', 'Sulfanilamide, CAM-B3LYP'),
                                ('sulfanilamide_pcm_d3/tddft_camb3lyp_eq', 'Sulfanilamide, CAM-B3LYP (eq.)')])]:
    rows = []
    for key_, flab in cases:
        t = D6['tddft'][key_]
        for s_ in t['states'][:10]:
            rows.append([flab if s_['state'] == 1 else '', f"S{s_['state']}", f(s_['E_eV'], 3), f(s_['nm'], 0), f(s_['f'], 4),
                         ', '.join(f"{a} ({100 * w:.0f} %)" for a, w in s_['dominant'][:2]).replace('->', '→'),
                         f(s_['CT'], 2) if s_['CT'] is not None else '–', f(s_['D_CT'], 2) if s_['D_CT'] is not None else '–'])
    table(ST / fn, ['Functional', 'State', '*E* (eV)', '*λ* (nm)', '*f*', 'Main configurations', 'CT number', '*D*~CT~ (Å)'],
          [':---', ':---', '---:', '---:', '---:', ':---', '---:', '---:'], rows,
          f'Lowest ten singlet excitations of the {lab} in IEF-PCM (water, non-equilibrium unless stated), TD-DFT/6-311++G(d,p) '
          'at the B3LYP-D3(BJ) geometry. H = HOMO, L = LUMO of the functional used. CT number: fraction of the excitation '
          'transferring charge between cation and anion; *D*~CT~: hole–electron centroid distance.', label,
          notes=[{'s-tddft-pair': 'All 20 roots converged; states 11–20 are in the result files.',
                  's-tddft-fragments': 'Twenty roots were requested and converged in each run; only the first ten are listed.'}[label]], doc='si', section='S11')

# S24 NLO details
rows = []
for sp, mod, k in NL:
    rows.append([sp, mod, f(v(f'nlo.{k}.mu'), 3), f(v(f'nlo.{k}.alpha'), 2), f(v(f'nlo.{k}.alpha_A3'), 2), f(v(f'nlo.{k}.dalpha'), 2),
                 f(v(f'nlo.{k}.beta'), 2), f(v(f'nlo.{k}.beta_B'), 2), f(v(f'nlo.{k}.beta_par'), 2), f(v(f'nlo.{k}.beta_esu'), 3)])
table(ST / 'tableS24_nlo.md', ['Species', 'Model', '*μ* (D)', '*ᾱ* (a.u.)', '*ᾱ* (Å³)', 'Δ*α* (a.u.)', '*β*~tot~ T (a.u.)', '*β*~tot~ B (a.u.)',
                               '*β*~∥~ (a.u.)', '*β*~tot~ (10⁻³⁰ esu, T)'],
      [':---', ':---', '---:', '---:', '---:', '---:', '---:', '---:', '---:', '---:'], rows,
      'Finite-field NLO properties, B3LYP-D3(BJ)/6-311++G(d,p): fields ±0.001 and ±0.002 a.u. along each axis and '
      'each pair of axes, one Richardson extrapolation; *β*~∥~ is the component along the dipole.', 's-nlo',
      notes=[f"The Richardson step changes *β*~tot~ of the IEF-PCM ion pair by {f(v('nlo.richardson_effect'), 1)} a.u.; the "
             'finite-field polarisability agrees with an analytic CPKS calculation to 2 × 10⁻⁷ (relative).'], doc='si', section='S12')

# S25 docking validation (all runs)
rows = []
for job, s_ in vs.items():
    st = s_['stats']
    g = J(DOCK / 'analysis/gate_G2.json')['criteria'].get(job, {})
    bv = lambda q: F[f"bulk.dockval.{job.replace('/', '.')}.{q}"]['value'] if st.get(q) is not None else None   # full precision
    rows.append([job.replace('__', ' '), f(bv('rmsd_top_mean'), 2, exact=True), f(bv('rmsd_top_sd'), 2), f(st.get('success_fraction_1A'), 1),
                 f(st.get('success_fraction_2A'), 1), f(bv('score_mean'), 2, exact=True), f(st.get('posebusters_physics_pass_fraction'), 1),
                 ('pass' if g.get('passes_G2') else 'fail') if g else '–', ('pass' if g.get('passes_2A_criterion') else 'fail') if g else '–'])
table(ST / 'tableS25_dock_validation.md', ['Run', 'RMSD (Å)', 'SD (Å)', 'f ≤ 1 Å', 'f ≤ 2 Å', 'Score (kcal mol⁻¹)', 'PoseBusters', 'Gate 1.0 Å', '2.0 Å'],
      [':---', '---:', '---:', '---:', '---:', '---:', '---:', ':---', ':---'], rows,
      'All validation and sensitivity docking runs: top-pose RMSD to the crystal ligand (mean, '
      'SD over seeds), success fractions, mean score and PoseBusters physics pass fraction.', 's-dock-validation',
      notes=['Gate 1.0 Å and 2.0 Å are the strict gate and the 2.0 Å criterion of this section. For the hCA II structures a run '
             'passes only if the zinc and Thr199 contacts of the site check are met as well, so a run with a fraction of 1.0 '
             'within 2 Å can fail the 2.0 Å criterion (Vina and Vinardo in hCA II; 6RL9 with Vina: Zn–N '
             f"{f(v('dock.val.T3_6RL9.sfa_neutral__vina.ZnN'), 2)} Å against {f(v('dock.gate.ZnN_max'), 1)} Å). The suffix of a run name gives "
             'the setting that the run changes (exhaustiveness, number of seeds, box edge, start conformer, loop build, '
             'His64 tautomer, sulfate removed, repeated run), as described in the text.'], doc='si', section='S13')

# S26 all production scores
rows = []
for key_, ro in ps['rank_order'].items():
    t, sc = key_.split('/')
    for rr in ro['ranking']:
        bk = f"bulk.dockprod.{t}.{sc}.{rr['species']}"
        rows.append([t.split('_')[1], sc, SPN[rr['species']], f(v(bk + '.score_mean'), 3, exact=True), f(v(bk + '.score_sd'), 3), f(rr['best'], 3), rr['rank'],
                     f(v(bk + '.LE'), 3), f(v(bk + '.FQ'), 2)])
# Units in the header cells, not in the caption: a single caption unit is inherited by every numeric column in
# numbers_audit.py, which then tried to read the dimensionless LE and FQ columns as kcal mol⁻¹ (29 Sep 2026).
table(ST / 'tableS26_dock_scores.md', ['Receptor', 'Scoring', 'Ligand', 'Mean (kcal mol⁻¹)', 'SD (kcal mol⁻¹)', 'Best (kcal mol⁻¹)', 'Rank', 'LE', 'FQ'],
      [':---', ':---', ':---', '---:', '---:', '---:', '---:', '---:', '---:'], rows,
      'Production docking scores of every species against every validated receptor. Top pose of each seed, with ligand '
      'efficiency LE = −score/heavy atoms and the size-corrected fit quality FQ (Reynolds scale).', 's-dock-scores',
      notes=['LE and FQ are size-biased heuristics (r(score, heavy atoms) = −0.87 with Vina) and are not used to rank potency.'], doc='si', section='S13')

# S26b relaxed torsion scans (referee M9; facts section 35) -- printed as Table S30
TSJ = J(RES / 'torsion_scan/torsion_strain.json')
PTS = {}
for sp in ('sfa_neutral', 'sfa_anion'):
    for p_ in D13['m9'][sp]['points']:
        if p_['branch'] != 'C3':
            continue
        PTS.setdefault(round(p_['target'], 1), {})[sp] = p_


def _lab(labels):
    l_ = labels[0]
    if l_.startswith('grid'):
        return 'grid'
    if 'bound' in l_:
        return l_.split(' ')[0] + ' crystal ligand' + (' (sulfamethoxazole)' if 'sulfamethoxazole' in l_ else ' (sulfanilamide)')
    mm = re.match(r'^(\w{4}) \([^)]*\), (Vina|Vinardo|AD4Zn), (representative|mean)', l_)
    return f"docked, {mm.group(1)} {mm.group(2)} ({'mean of the top poses' if mm.group(3) == 'mean' else 'representative pose'})"


rows = [[f(v('tors.ref.sfa_neutral.tau'), 1), 'IEF-PCM minimum of sulfanilamide', '0', '–'],
        [f(v('tors.ref.sfa_anion.tau'), 1), f"IEF-PCM minimum of the sulfonamidate, re-optimised ({f(v('frag.sulfonamidate_pcm_d3.tau_fold'), 1)}° in the docking input)", '–', '0'],
        [f(v('xray.tau_C3C4S1N2'), 1), 'SAPH crystal (cation)^a^', f(v('tors.strain.sfa_neutral.saph_crystal'), 2), f(v('tors.strain.sfa_anion.saph_crystal'), 2)]]
for t_ in sorted(PTS):
    d_ = PTS[t_]
    labs = sorted({_lab(q['labels']) for q in d_.values()})
    kk = f"t{int(round(t_ * 10)):04d}"
    rows.append([f(t_, 1), '; '.join(labs)] + [f(v(f'tors.scan.{sp}.{kk}'), 2) if sp in d_ else '–' for sp in ('sfa_neutral', 'sfa_anion')])
rows = sorted(rows, key=lambda r_: float(r_[0]))
table(ST / 'tableS26b_torsion_scan.md', ['*τ*~fold~ (°)', 'Point', 'Sulfanilamide (kcal mol⁻¹)', 'Sulfonamidate (kcal mol⁻¹)'],
      ['---:', ':---', '---:', '---:'], rows,
      'Relaxed torsional strain of neutral sulfanilamide and the sulfonamidate: C3–C4–S1–N2 held, everything else optimised, '
      'B3LYP-D3(BJ)/6-311++G(d,p) in IEF-PCM (water); electronic energy relative to the unconstrained minimum of the same '
      'species. *τ*~fold~ = max(|C3–C4–S1–N2|, |C5–C4–S1–N2|).', 's-torsion-scan',
      notes=['^a^ Interpolated (cubic spline through the computed points and the minimum; the spline reproduces the exact off-grid '
             f"points within {f(v('tors.scan.sfa_neutral.interp_check'), 3)} kcal mol⁻¹); every other value is a constrained optimisation. "
             f"All {v('tors.scan.sfa_neutral.n_converged')} sulfanilamide and {v('tors.scan.sfa_anion.n_converged')} sulfonamidate points converged. A point on the "
             f"other ring-flip branch (C5–C4–S1–N2 held at 150°) differs from the C3 branch by {f(v('tors.scan.sfa_neutral.mirror_diff'), 3)} kcal mol⁻¹. "
             'The sulfonamidate profile is flat: its re-optimised minimum is not the lowest point of the scan, which lies at '
             f"{f(v('tors.scan.sfa_anion.tau_at_min'), 0)}° ({f(v('tors.scan.sfa_anion.min'), 2)} kcal mol⁻¹)."], doc='si', section='S13')

# S27 torsions and strain
# Numbers audit round 1 (29 Sep 2026): the column "vs optimum with the same N–H" (facts section 31, dock.strainref.*)
# was added because the rigid single points against the optimised minimum carry a common offset from the MMFF94s
# placement of the N–H hydrogens (3.0 kcal mol⁻¹ for sulfanilamide, 4.5 for the sulfonamidate; no reference single
# point for the cation). The old header "N–H re-optimised" was wrong: the hydrogens were placed by MMFF94s.
def _sref(k_):
    return f(F[k_]['value'], 2) if F.get(k_, {}).get('value') is not None else '–'


# 29 Sep 2026 (VM burst 3): the rigid DFT single points (facts sections 29 and 31) are retired; the DFT column is now the
# relaxed torsional strain of the scan (facts section 35) at the torsion of every docked seed, and the crystal-ligand
# torsion of each receptor is set beside the docked ones (referee chem M9c, dock M6).
XTAL = {'3TZF': (v('dock.cocrystal.3TZF.tau'), ''), '6RL9': (v('dock.cocrystal.6RL9.tau'), ''), '3HS4': (v('dock.cocrystal.6RL9.tau'), '^b^'),
        '7S2K': (None, '^c^')}
rows = [['Crystal (SAPH, CIF)', '–', '–', f(v('xray.tau_C3C4S1N2'), 1), '–', f"{f(v('tors.strain.sfa_neutral.saph_crystal'), 2)} / {f(v('tors.strain.sfa_anion.saph_crystal'), 2)}^a^", '–']]
for sp, c in [('sfa_neutral', 'Sulfanilamide'), ('sfa_anion', 'Sulfonamidate'), ('sfa_cation', 'Cation')]:
    rows.append([f'DFT, {c.lower()} (IEF-PCM)', '–', '–', f(v(f'dock.tors.dft.{sp}'), 1), '–', '0 (minimum)' if sp != 'sfa_cation' else '–', '–'])
for job, j in ps['jobs'].items():
    if j['ligand'].startswith('sfa_') and j.get('torsion') and not j.get('diagnostic'):
        tt = j['torsion']
        tg_ = job.split('/')[0].split('_')[1]
        kx = f"tors.seeds.{j['ligand']}.{tg_}.{j['scoring']}"
        xt, xm = XTAL[tg_]
        # closing items (29 Sep 2026): tau_fold and the DFT strain from the raw pose coordinates (facts section 37,
        # dock.tors.raw.*), which the text and Fig. 10 use; the rebuilt-pose values of production_summary.json differ by
        # up to 0.2 deg and give the same strain ranges at two decimals.
        kr = f"dock.tors.raw.{tg_}.{j['ligand']}.{j['scoring']}"
        tr_ = v(kr)
        ks_ = kr + '.strain'
        rows.append([tg_ + ', ' + {'sfa_neutral': 'sulfanilamide', 'sfa_anion': 'sulfonamidate', 'sfa_cation': 'cation'}[j['ligand']],
                     j['scoring'], j['n_seeds'], _rng({'min': tr_[0], 'max': tr_[1]}, 1),  # prose audit burst 3: a single value when min = max
                     (f(xt, 1) if xt else '–') + xm,
                     ((f(v(ks_)[0], 2) if f(v(ks_)[0], 2) == f(v(ks_)[1], 2) else f"{f(v(ks_)[0], 2)} to {f(v(ks_)[1], 2)}") if ks_ in F else '–'),
                     _rng(tt['strain_mmff_relaxed'], 2)])
table(ST / 'tableS27_dock_torsion.md', ['Structure', 'Scoring', 'Seeds', '*τ*~fold~ of the top poses (°)', 'Crystal ligand *τ*~fold~ (°)',
                                        'Relaxed DFT strain (kcal mol⁻¹)', 'MMFF94s strain, relaxed (kcal mol⁻¹)'],
      [':---------------------------------', ':----------------', '------------:', '----------------------------------:', '------------------------------:', '----------------------------------:', '----------------------------------:'], rows,
      'The sulfonamide torsion from the crystal to the docked poses, with its torsional strain. *τ*~fold~ = '
      'max(|C3–C4–S1–N2|, |C5–C4–S1–N2|) in the SAPH crystal, the DFT optima and the docked top poses (range over seeds), beside the torsion of the ligand in the deposited complex of each receptor, with the '
      'relaxed B3LYP-D3(BJ)/6-311++G(d,p) IEF-PCM (water) torsional strain at the torsions of the docked poses (range over seeds) '
      'and the MMFF94s estimate.', 's-dock-torsion',
      notes=['*τ*~fold~ of the docked poses is measured on the top pose of each seed as written by the docking program. '
             'Relaxed DFT strain: the profile of @tbl:s-torsion-scan (C3–C4–S1–N2 held, everything else optimised) evaluated at the '
             '*τ*~fold~ of each seed; not computed for the cation. ^a^ Sulfanilamide / sulfonamidate at the crystal torsion. '
             f"^b^ The 6RL9 sulfanilamide, carried into 3HS4 by superposition, is the cross-docking reference. ^c^ The 7S2K ligand, "
             'dihydropteroate, has no sulfonamide. MMFF94s: minimised with the torsion restrained, relative to the global MMFF94s '
             f"minimum; it follows the DFT strain of sulfanilamide within {f(max(abs(x_) for x_ in v('tors.mmff_minus_dft.sfa_neutral')), 2)} kcal mol⁻¹ "
             f"at the representative poses but overstates that of the sulfonamidate by {f(v('tors.mmff_minus_dft.sfa_anion')[0], 2)}–"
             f"{f(v('tors.mmff_minus_dft.sfa_anion')[1], 2)} kcal mol⁻¹."],
      doc='si', section='S13')

# S28 contacts
pcj = J(DOCK / 'results/production_contacts.json')
rows = []
for job in ['T5_7S2K/sfa_neutral__vina__prod', 'T5_7S2K/sfa_anion__vina__prod', 'T1_3TZF/sfa_neutral__vina__prod', 'T1_3TZF/sfa_anion__vina__prod',
            'T3_6RL9/sfa_neutral__ad4zn__prod', 'T4_3HS4/sfa_neutral__ad4zn__prod', 'T3_6RL9/sfa_neutral__vina__prod', 'T4_3HS4/sfa_neutral__vina__prod']:
    c1 = pcj[job]['1']
    rn = lambda k: {'HH2': 'pterin (HH2)', 'POP': 'PPi', 'HOH': 'water'}.get(k[:3], k[:1] + k[1:3].lower() + k[3:]) if not k[:3].isdigit() else k
    pol = '; '.join(f"{rn(p_['residue'])} {p_['receptor_atom']}···{p_['ligand_atom'].split('#')[0]} {p_['distance']:.2f}" for p_ in c1['polar_pairs_NO_within_3.5A'][:4])
    zn = c1.get('zinc') or {}
    rows.append([job.split('/')[0].split('_')[1], job.split('/')[1].replace('__prod', '').replace('__', ', '), f"{c1['population']}", pol or '–',
                 f"{zn['nearest_ligand_atom'].split('#')[0]} {zn['distance']:.2f}" if zn.get('nearest_ligand_atom') else '–'])
table(ST / 'tableS28_dock_contacts.md', ['Receptor', 'Species, scoring', 'Cluster population', 'Polar contacts ≤ 3.5 Å (Å)', 'Nearest atom to Zn (Å)'],
      [':---', ':---', '---:', ':---', ':---'], rows,
      'Contacts of the representative top pose of the largest cluster for the sulfanilamide species (heavy-atom distances; '
      'PLIP/ProLIF/geometric).', 's-dock-contacts', doc='si', section='S13')

# ------------------------------------------------------------------------------------------------ control receptor (4 Oct 2026)
# S29/S30 by file name, printed as Tables S33 and S34 (after the contacts table): the ATP site of E. coli DNA gyrase B
# (4DUH), not a target of sulfanilamide, docked with the protocol of the five receptors (facts section 38; dock/control_gyrb.py).
CSJ = J(DOCK / 'control/control_summary.json')


def _cv(key, d=2, sd=None):
    """A control fact at d decimals, with '± SD' when the SD key is given; '–' when the job was not run."""
    if key not in F:
        return '–'
    return f(v(key), d, exact=True) + (f' ± {f(v(sd), 2)}' if sd and sd in F else '')


def _crange(key, d=2):
    """'min–max' of a two-member control fact. The members are exact (integers, or differences of the three-decimal
    scores that Vina prints), so a value ending in 5 is a true tie and is rounded half to even without a warning."""
    if key not in F:
        return '–'
    a_, b_ = f(v(key)[0], d, exact=True), f(v(key)[1], d, exact=True)
    return a_ if a_ == b_ else f'{a_}–{b_}'


rows = []
for lab_, tag_, sc_ in [('As prepared^a^', 'primary', 'vina'), ('As prepared^a^', 'primary', 'vinardo'),
                        ('Residues 109–116 built from a second seed', 'second', 'vina'), ('Built residues 109–116 removed', 'noloop', 'vina'),
                        ('Ordered crystal water kept^b^', 'water', 'vina')]:
    b_ = f'ctrl.val.{tag_}.{sc_}'
    if b_ + '.rmsd' not in F:
        continue
    kk_ = {'primary': f'ctrl.score.{{}}.{sc_}', 'water': 'ctrl.wat.{}.vina', 'second': 'ctrl.var.second.{}.score', 'noloop': 'ctrl.var.noloop.{}.score'}[tag_]
    ks_ = {'primary': f'ctrl.score.{{}}.{sc_}.sd', 'water': 'ctrl.wat.{}.vina.sd', 'second': 'ctrl.var.second.{}.sd', 'noloop': 'ctrl.var.noloop.{}.sd'}[tag_]
    ki_ = {'primary': f'ctrl.score.sfa_neutral.{sc_}.in_site', 'water': 'ctrl.wat.sfa_neutral.vina.in_site', 'second': 'ctrl.var.second.sfa_neutral.in_site',
           'noloop': 'ctrl.var.noloop.sfa_neutral.in_site'}[tag_]
    nn_ = v(b_ + '.nn.n_seeds')
    rows.append([lab_, {'vina': 'Vina', 'vinardo': 'Vinardo'}[sc_], f"{f(v(b_ + '.rmsd'), 2, exact=True)} ± {f(v(b_ + '.rmsd_sd'), 2)}",
                 f"{f(v(b_ + '.f1'), 1)} ({f(v(b_ + '.f2'), 1)})", _cv(b_ + '.score', 2, b_ + '.score_sd'),
                 (_crange(b_ + '.nn.rank', 0) if nn_ else '–'), (_crange(b_ + '.nn.gap', 1) if nn_ else '–'), _crange(b_ + '.top_waters', 0),
                 _cv(kk_.format('sfa_neutral'), 2, ks_.format('sfa_neutral')), (f"{v(ki_)}/{v('protocol.dock_seeds')}" if ki_ in F else '–'),
                 _cv(kk_.format('cfx_zwit'), 2, ks_.format('cfx_zwit'))])
table(ST / 'tableS29_control_receptors.md',
      ['Receptor', 'Scoring', 'Inhibitor RMSD (Å)', '*f* ≤ 1 Å (≤ 2 Å)', 'Inhibitor score (kcal mol⁻¹)', 'Rank of the near-native mode',
       'Its score behind the top pose (kcal mol⁻¹)', 'Waters under the top pose', 'Sulfanilamide (kcal mol⁻¹)', 'Its seeds in the cavity^c^',
       'Ciprofloxacin (kcal mol⁻¹)'],
      [':----------------', ':-------', '-------------:', '---------:', '--------------:', '-------:', '-----------:', '-------:', '--------------:', '-------:', '--------------:'], rows,
      'Control receptor, the ATP site of *Escherichia coli* DNA gyrase B (4DUH, chain A): redocking of the crystal inhibitor in every '
      'variant of the receptor, with the scores of sulfanilamide and ciprofloxacin in the same variant. RMSD is the symmetry-corrected heavy-atom RMSD of '
      'the top pose from the crystal pose in the receptor frame, mean ± SD over 10 seeds, with the fraction of seeds within 1 Å (2 Å); the '
      'near-native mode is the first mode of a seed within 1.0 Å of the crystal pose, given by its rank and by its score minus the top score '
      '(range over the seeds); "waters under the top pose" counts the crystal waters of chain A whose oxygen lies within '
      f"{f(v('ctrl.water_clash'), 1)} Å of a heavy atom of the top pose. Scores are means ± SD over the seeds.", 's-control-receptors',
      notes=['^a^ Prepared like the five receptors of Section 2.5 of the main text, with the unobserved lid residues Tyr109–His116 built by PDBFixer. '
             f"^b^ The receptor as prepared with {v('ctrl.water.n_kept')} of the {v('ctrl.water.n_chain')} waters of chain A added as fixed oxygen atoms "
             f"({v('ctrl.water.n_in_box')} inside the box); a diagnostic, since the protocol keeps no such water. "
             f"^c^ Seeds whose top pose of sulfanilamide has its centroid within {f(v('ctrl.site_cut'), 1)} Å of a heavy atom of the crystal inhibitor; its "
             'other top poses lie under the ATP lid, on four crystal waters. '
             'Every top pose of this table passes the PoseBusters physics checks.'],
      doc='si', section='S13')

# note c of the table above: every sulfanilamide top pose outside the cavity of the inhibitor overlaps four crystal waters
for k_, q_ in CSJ['jobs'].items():
    if '/sfa_neutral__' in k_:
        assert all(r_['crystal_waters_overlapped'] == 4 for r_ in q_['site']['per_seed'] if not r_['in_site']), k_
CTN = dict(SPN, rli_anion='4DUH inhibitor (crystal ligand, RLI)', cfx_zwit='Ciprofloxacin')
rows = []
for r_ in CSJ['score_table']['vina']['ranking']:
    sp_ = r_['species']
    b_, w_ = f'ctrl.score.{sp_}.vina', f'ctrl.wat.{sp_}.vina'
    tg_ = [(f(v(f'bulk.dockprod.{t_}.vina.{sp_}.score_mean'), 2, exact=True) if f'bulk.dockprod.{t_}.vina.{sp_}.score_mean' in F else '–')
           for t_ in ('T1_3TZF', 'T5_7S2K', 'T3_6RL9', 'T4_3HS4')]
    rows.append([CTN[sp_], v(b_ + '.heavy'), _cv(b_, 2, b_ + '.sd'), _crange(b_ + '.waters_range', 0), f"{v(b_ + '.in_site')}/{v(b_ + '.n_seeds')}",
                 _cv(w_, 2, w_ + '.sd'), (f"{v(w_ + '.in_site')}/{v(w_ + '.n_seeds')}" if w_ + '.in_site' in F else '–')] + tg_)
_dn = lambda k_: f"{f(v(k_ + '.mean'), 2, sign=True)} ({f(v(k_ + '.range')[0], 2, sign=True)} to {f(v(k_ + '.range')[1], 2, sign=True)})"  # noqa: E731
_z8 = 'ctrl.size8.vina'
table(ST / 'tableS30_control_scores.md',
      ['Ligand', 'Heavy atoms', '4DUH as prepared (kcal mol⁻¹)', 'Waters under the pose^a^', 'In the cavity^b^', '4DUH with crystal water (kcal mol⁻¹)', 'In the cavity, water^b^',
       '3TZF (kcal mol⁻¹)', '7S2K (kcal mol⁻¹)', '6RL9 (kcal mol⁻¹)', '3HS4 (kcal mol⁻¹)'],
      [':----------------------------', '------:', '--------------:', '-------:', '------:', '--------------:', '------:', '--------:', '--------:', '--------:', '--------:'], rows,
      'Vina scores in the control receptor (4DUH) beside the scores of the same species in the four production structures. Mean ± SD of the top-pose '
      'score over the seeds in the receptor as prepared and in the receptor with the ordered crystal water kept, ordered by the first; the scores in '
      '3TZF, 7S2K, 6RL9 and 3HS4 are the production means of @tbl:s-dock-scores with Vina.', 's-control-scores',
      notes=[f"^a^ Crystal waters of chain A whose oxygen lies within {f(v('ctrl.water_clash'), 1)} Å of a heavy atom of the top pose in the receptor as prepared (range over the seeds). "
             f"^b^ Seeds whose top pose has its centroid within {f(v('ctrl.site_cut'), 1)} Å of a heavy atom of the crystal inhibitor. The count says where a pose lies, not that it is "
             f"correct: the top pose of the inhibitor as prepared counts in all seeds, {f(v('ctrl.val.primary.vina.rmsd'), 2, exact=True)} Å from its crystal pose. "
             + (f"Over the {v(_z8 + '.n')} species docked into every receptor, the Pearson correlation between the score and the number of heavy atoms is "
                f"{f(v(_z8 + '.4DUH.r'), 2)} in 4DUH as prepared"
                + (f" ({f(v('ctrl.size8.vina_water.4DUH.r'), 2)} with the crystal water)" if 'ctrl.size8.vina_water.4DUH.r' in F else '')
                + f", {f(v(_z8 + '.6RL9.r'), 2)} in 6RL9, {f(v(_z8 + '.3HS4.r'), 2)} in 3HS4, {f(v(_z8 + '.3TZF.r'), 2)} in 3TZF and {f(v(_z8 + '.7S2K.r'), 2)} in 7S2K; "
                f"over the {v('ctrl.size.vina.n')} species in 4DUH as prepared it is {f(v('ctrl.size.vina.r'), 2)}. " if _z8 + '.n' in F else '')
             + (f"Mean difference of the score in 4DUH as prepared from that in each production structure over these {v('ctrl.cmp.vina.3TZF.n')} species, in kcal mol⁻¹ "
                f"(range in parentheses): 3TZF {_dn('ctrl.cmp.vina.3TZF')}, 7S2K {_dn('ctrl.cmp.vina.7S2K')}, 6RL9 {_dn('ctrl.cmp.vina.6RL9')}, 3HS4 {_dn('ctrl.cmp.vina.3HS4')}"
                if 'ctrl.cmp.vina.3TZF.mean' in F else '')
             + (f"; with the crystal water, {_dn('ctrl.cmp.vina_water.3TZF')}, {_dn('ctrl.cmp.vina_water.7S2K')}, {_dn('ctrl.cmp.vina_water.6RL9')} and {_dn('ctrl.cmp.vina_water.3HS4')}. "
                if 'ctrl.cmp.vina_water.3TZF.mean' in F else ('. ' if 'ctrl.cmp.vina.3TZF.mean' in F else ''))
             + ('A positive difference is a weaker score in the control receptor.' if 'ctrl.cmp.vina.3TZF.mean' in F else '')], doc='si', section='S13')

# ------------------------------------------------------------------------------------------------ SI completion (29 Sep 2026)
# S9b condensed Fukui functions (facts section 24; paper/si/compute/fukui.py)
def fk_rows(sp, label):
    fk = J(RES / {'sfa': 'sulfanilamide_pcm_d3', 'cation': 'cation_pcm_d3'}[sp] / 'fukui/fukui.json')
    out = []
    for i, a_ in enumerate(fk['atoms']):
        cells = [f(v(f'fukui.{sp}.{s_}.{k_}.{a_}'), 3) for s_ in ('npa', 'hirsh') for k_ in ('fp', 'fm', 'f0')]
        out.append([label if i == 0 else '', a_] + cells)
    for g, gl in [('amino', 'Σ NH₂' if sp == 'sfa' else 'Σ NH₃⁺'), ('ring', 'Σ C₆H₄'), ('sulfonamide', 'Σ SO₂NH₂')]:
        out.append(['', gl] + [f(v(f'fukui.{sp}.{s_}.{k_}.grp.{g}'), 2) for s_ in ('npa', 'hirsh') for k_ in ('fp', 'fm', 'f0')])
    return out


if (RES / 'cation_pcm_d3/fukui/fukui.json').exists() and (RES / 'sulfanilamide_pcm_d3/fukui/fukui.json').exists():
    rows = fk_rows('sfa', 'Sulfanilamide') + fk_rows('cation', 'Cation')
    table(ST / 'tableS_fukui.md', ['Species', 'Atom', '*f*⁺ (NPA)', '*f*⁻ (NPA)', '*f*⁰ (NPA)', '*f*⁺ (Hirshfeld)', '*f*⁻ (Hirshfeld)', '*f*⁰ (Hirshfeld)'],
          [':---', ':---', '---:', '---:', '---:', '---:', '---:', '---:'], rows,
          'Condensed Fukui functions of neutral sulfanilamide and the 4-sulfamoylanilinium cation in IEF-PCM (water), '
          'B3LYP-D3(BJ)/6-311++G(d,p), each at its own IEF-PCM optimum. Finite differences of atomic charges *q* of the '
          '*N*, *N* + 1 and *N* − 1 electron systems at fixed geometry (restricted Kohn–Sham for *N*, unrestricted doublets for '
          '*N* ± 1, each state solvated at equilibrium): *f*⁺ = *q*(*N*) − *q*(*N* + 1), *f*⁻ = *q*(*N* − 1) − *q*(*N*), '
          '*f*⁰ = [*q*(*N* − 1) − *q*(*N* + 1)]/2 [@parr1984density; @yang1986use]. NPA charges from PySCF natural atomic '
          'orbitals of each state [@reed1985natural]; Hirshfeld charges [@hirshfeld1977bonded] with a promolecule of '
          'spherically averaged neutral atoms computed at the same level. Each column sums to 1 over the atoms of a species; '
          'Σ rows are group sums.', 's-fukui',
          notes=[f"The three total energies of each species equal the ΔSCF runs behind @tbl:s-descriptors (vertical *I* "
                 f"{f(v('fukui.sfa.pcm.I'), 3)} and *A* {f(v('fukui.sfa.pcm.A'), 3)} eV for sulfanilamide, {f(v('fukui.cation.pcm.I'), 3)} and "
                 f"{f(v('fukui.cation.pcm.A'), 3)} eV for the cation). ⟨*S*²⟩ of the doublets: {f(v('fukui.sfa.pcm.S2_Nm1'), 3)} and "
                 f"{f(v('fukui.sfa.pcm.S2_Np1'), 3)} (sulfanilamide *N* − 1, *N* + 1), {f(v('fukui.cation.pcm.S2_Nm1'), 3)} and "
                 f"{f(v('fukui.cation.pcm.S2_Np1'), 3)} (cation). Small negative values are an artefact of condensing finite "
                 "differences onto atoms."], doc='si', section='S9')

# S14 ADMET and drug-likeness (facts section 25; paper/si/compute/admet.py)
if (RES / 'sulfanilamide_pcm_d3/admet/admet.json').exists():
    def tr(ep):
        return '†' if v(f'admet.sfa.train.{ep}') else ''
    rows = [
        ['*Drug-likeness (RDKit)*', '', '', '', ''],
        ['Molecular weight (g mol⁻¹)', f(v('admet.sfa.MW'), 2), f(v('admet.cation.MW'), 2), '–', ''],
        ['cLogP (Wildman–Crippen)', f(v('admet.sfa.cLogP'), 2), f(v('admet.cation.cLogP'), 2), f"log *P* {f(v('admet.exp.logP'), 2)}", '[@knox2024drugbank]'],
        ['Molar refractivity (cm³ mol⁻¹)', f(v('admet.sfa.MR'), 1), f(v('admet.cation.MR'), 1), '–', ''],
        ['TPSA, N and O (Å²)', f(v('admet.sfa.TPSA'), 1), f(v('admet.cation.TPSA'), 1), '–', ''],
        ['H-bond donors / acceptors (NH + OH / N + O)', f"{v('admet.sfa.HBD')} / {v('admet.sfa.HBA')}", f"{v('admet.cation.HBD')} / {v('admet.cation.HBA')}", '–', ''],
        ['Rotatable bonds', v('admet.sfa.rotb'), v('admet.cation.rotb'), '–', ''],
        ['Atoms, including H', v('admet.sfa.atoms'), v('admet.cation.atoms'), '–', ''],
        ['Violations: Lipinski / Veber / Ghose / Egan', ' / '.join(str(v(f'admet.sfa.{r_}.n_viol')) for r_ in ('lipinski', 'veber', 'ghose', 'egan')),
         ' / '.join(str(v(f'admet.cation.{r_}.n_viol')) for r_ in ('lipinski', 'veber', 'ghose', 'egan')), '–', ''],
        ['*Absorption and distribution (ADMET-AI)*', '', '', '', ''],
        ['Aqueous solubility, log *S* (mol L⁻¹)', f(v('admet.sfa.ai.Solubility_AqSolDB'), 2) + tr('Solubility_AqSolDB'), f(v('admet.cation.ai.Solubility_AqSolDB'), 2),
         f"{f(v('admet.exp.logS'), 2)} ({f(v('admet.exp.solubility'), 1)} g L⁻¹ at {v('admet.exp.T_sol')} °C)", '[@knox2024drugbank; @sorkun2019aqsoldb]'],
        ['Lipophilicity, log *D*~7.4~', f(v('admet.sfa.ai.Lipophilicity_AstraZeneca'), 2), f(v('admet.cation.ai.Lipophilicity_AstraZeneca'), 2),
         f"log *P* {f(v('admet.exp.logP'), 2)}", '[@knox2024drugbank]'],
        ['Human intestinal absorption (probability)', f(v('admet.sfa.ai.HIA_Hou'), 2), f(v('admet.cation.ai.HIA_Hou'), 2),
         f"{v('admet.exp.oral_absorbed')[0]}–{v('admet.exp.oral_absorbed')[1]} % of an oral dose absorbed (sulfonamides)", '[@kim2025pubchem]'],
        ['Oral bioavailability (probability)', f(v('admet.sfa.ai.Bioavailability_Ma'), 2), f(v('admet.cation.ai.Bioavailability_Ma'), 2),
         f"blood level {v('admet.exp.blood_6h')[0]}–{v('admet.exp.blood_6h')[1]} mg per 100 mL {v('admet.exp.blood_6h_time')} h after an oral dose (no bioavailability measured)", '[@taylor1940comparative]'],
        ['Plasma protein binding (%)', f(v('admet.sfa.ai.PPBR_AZ'), 0), f(v('admet.cation.ai.PPBR_AZ'), 0),
         f"about {f(v('admet.exp.ppb.bound'), 0)} (human plasma, dialysis)", '[@davis1943binding]'],
        ['Blood–brain barrier penetration (probability)', f(v('admet.sfa.ai.BBB_Martins'), 2), f(v('admet.cation.ai.BBB_Martins'), 2), '–', ''],
        ['*Excretion and metabolism (ADMET-AI)*', '', '', '', ''],
        ['Half-life (h)', f(v('admet.sfa.ai.Half_Life_Obach'), 1), f(v('admet.cation.ai.Half_Life_Obach'), 1), f"≈{f(v('admet.exp.half_life'), 0)}", '[@koizumi1964absorption; @nogami1968studies]'],
        ['CYP2C9 / CYP3A4 inhibition (probability)', f"{f(v('admet.sfa.ai.CYP2C9_Veith'), 2)}{tr('CYP2C9_Veith')} / {f(v('admet.sfa.ai.CYP3A4_Veith'), 2)}{tr('CYP3A4_Veith')}",
         f"{f(v('admet.cation.ai.CYP2C9_Veith'), 2)} / {f(v('admet.cation.ai.CYP3A4_Veith'), 2)}", '–', ''],
        ['*Toxicity (ADMET-AI)*', '', '', '', ''],
        ['hERG blocking (probability)', f(v('admet.sfa.ai.hERG'), 2), f(v('admet.cation.ai.hERG'), 2), '–', ''],
        ['Ames mutagenicity (probability)', f(v('admet.sfa.ai.AMES'), 2) + tr('AMES'), f(v('admet.cation.ai.AMES'), 2), 'negative', '[@hansen2009benchmark]'],
        ['Drug-induced liver injury (probability)', f(v('admet.sfa.ai.DILI'), 2) + tr('DILI'), f(v('admet.cation.ai.DILI'), 2), 'positive (DILIst)', '[@thakkar2020drug]'],
        ['Rat oral LD~50~, −log(LD~50~/mol kg⁻¹)', f(v('admet.sfa.ai.LD50_Zhu'), 2) + tr('LD50_Zhu'), f(v('admet.cation.ai.LD50_Zhu'), 2),
         f"{f(v('admet.exp.LD50_rat_oral_log'), 3)} ({f(v('admet.exp.LD50_rat_oral'), 1)} g kg⁻¹)", '[@zhu2009quantitative]'],
        ['*Ionisation (experimental only)*', '', '', '', ''],
        ['pK~a~, sulfonamide N–H', '–', '–', f"{f(v('admet.exp.pKa_SO2NH2'), 2)}; {f(v('admet.exp.pKa_SO2NH2_20C'), 1)} ({v('admet.exp.T_pKa')} °C)", '[@kim2025pubchem; @knox2024drugbank]'],
        ['pK~a~, anilinium', '–', '–', f"{f(v('admet.exp.pKa_anilinium'), 2)} ({v('admet.exp.T_pKa_anilinium')} °C); {f(v('admet.exp.pKa_anilinium_klotz'), 2)} ({v('admet.exp.T_pKa_anilinium_klotz')} °C)", '[@greenbaum1955comparison; @klotz1945isoelectric; @zheng2026digitized]'],
    ]
    table(ST / 'tableS_admet.md', ['Property', 'Sulfanilamide', 'Cation (reference)', 'Experimental, sulfanilamide', 'Source'],
          [':---', '---:', '---:', ':---', ':---'], rows,
          'Drug-likeness descriptors and ADMET-AI predictions for neutral sulfanilamide, the species present at pH 7.4, with '
          'the 4-sulfamoylanilinium cation for reference, against experimental values for sulfanilamide. Descriptors from RDKit '
          '(cLogP and molar refractivity [@wildman1999prediction], TPSA over N and O [@ertl2000fast]); filters of Lipinski '
          '[@lipinski2001experimental], Veber [@veber2002molecular], Ghose [@ghose1999knowledge] and Egan [@egan2000prediction]; '
          'ADMET-AI [@swanson2024admet] probabilities are classifier outputs, not measured frequencies. † Sulfanilamide is in the '
          'training data of this endpoint, so the prediction is not independent of the experimental value. The cation lies '
          'outside the training domain of ADMET-AI.', 's-admet',
          notes=['Sulfanilamide fails only the Ghose atom count (19 atoms, range 20–70); the cation fails only the Ghose logP range. '
                 'The oral-absorption range is a class statement for sulfonamides; the blood level was measured 6 h after an oral dose. '
                 'Plasma protein binding: bound/free ratio 4.6 % per g/100 mL of protein at 10 mg/100 mL drug, evaluated at 6 g/100 mL '
                 'protein. Half-life: ln 2/*k* from the overall elimination rate constant of sulfanilamide in man (about 0.067 h⁻¹, urinary excretion after a '
                 '1 g oral dose); the conversion is ours, neither source prints a half-life. Anilinium p*K*~a~: UV absorbance in dilute acid, '
                 'as evaluated in the IUPAC compilation.'], doc='si', section='S14')

# ------------------------------------------------------------------------------------------------ docking tests and MD (8 Oct 2026)
# S31 to S35 by file name, printed as Tables S36 to S40 (after the ADMET table): the two docking tests in the hCA II
# receptors (Section S15: Glu106 as the acid of the printed results or as a carboxylate; four crystal waters of 6RL9
# kept) and the molecular dynamics of the hCA II-sulfanilamide complex (Section S16: settings, the eleven runs, the
# tests of the force field).  Every cell is a fact of section 40 (facts_md.py; derive/d17_md.py).
def _n(x_):
    """An integer with a thousands separator."""
    return f'{int(x_):,}'


def _pm(k_, d_=2):
    return f"{f(v(k_ + '.rmsd'), d_)} ± {f(v(k_ + '.rmsd_sd'), d_)}"


def _sc(k_, d_=2):
    return f"{f(v(k_ + '.score'), d_, exact=True)} ± {f(v(k_ + '.score_sd'), d_)}"


def _gates(k_):
    return f"{v(k_ + '.gate')} / {v(k_ + '.crit2A')}" if k_ + '.gate' in F else '–'


_LG = {'sfa_neutral': 'Sulfanilamide', 'sfa_anion': 'Sulfonamidate', 'azm_neutral': 'Acetazolamide', 'azm_anion': 'Acetazolamide anion'}
_SCN = {'ad4zn': 'AD4Zn', 'vina': 'Vina', 'vinardo': 'Vinardo'}
_GLU = {'acid': 'acid', 'carb': 'carboxylate'}


def _gate2(k_):
    return [v(k_ + '.gate'), v(k_ + '.crit2A')] if k_ + '.gate' in F else ['–', '–']


rows, _sds = [], {'ad4zn': [], 'vina': [], 'vinardo': []}
for lig_, sc_, kind_ in (('sfa_neutral', 'ad4zn', 'prod'), ('sfa_neutral', 'ad4zn', 'val'), ('sfa_anion', 'ad4zn', 'val'),
                         ('sfa_neutral', 'vina', 'val'), ('sfa_neutral', 'vinardo', 'val')):
    first_ = True
    recs_ = [('acid', 'no', f'dock.glu106.6RL9.{lig_}.{sc_}.{kind_}.acid', f'dock.wat4.dry.acid.{lig_}.{sc_}.{kind_}'),
             ('carb', 'no', f'dock.glu106.6RL9.{lig_}.{sc_}.{kind_}.carb', f'dock.wat4.dry.carb.{lig_}.{sc_}.{kind_}'),
             ('acid', 'yes', f'dock.wat4.acid.{lig_}.{sc_}.{kind_}', f'dock.wat4.acid.{lig_}.{sc_}.{kind_}'),
             ('carb', 'yes', f'dock.wat4.carb.{lig_}.{sc_}.{kind_}', f'dock.wat4.carb.{lig_}.{sc_}.{kind_}')]
    if (lig_, sc_, kind_) == ('sfa_neutral', 'ad4zn', 'prod'):
        recs_ += [('acid', 'yes^a^', 'dock.wat4.rev.sfa_neutral.ad4zn.prod', 'dock.wat4.rev.sfa_neutral.ad4zn.prod'),
                  ('acid', 'yes^b^', 'dock.wat4.bare.sfa_neutral.ad4zn.prod', 'dock.wat4.bare.sfa_neutral.ad4zn.prod')]
    for glu_, wat_, k_, km_ in recs_:
        lab_ = f"{_LG[lig_]}, {_SCN[sc_]}, {'DFT geometry' if kind_ == 'prod' else 'ETKDG'}, {v(k_ + '.n_seeds')} seeds" if first_ else ''
        first_ = False
        _sds[sc_].append(v(k_ + '.score_sd'))
        rows.append([lab_, _GLU[glu_], wat_, _pm(k_), f"{v(k_ + '.n_le_1')} of {v(k_ + '.n_seeds')}", v(km_ + '.modes_rot'), f(v(k_ + '.score'), 3, exact=True),
                     f(v(k_ + '.zn_n'), 2)] + _gate2(k_))
table(ST / 'tableS31_dock_glu106_water.md',
      ['Ligand, scoring, input', 'Glu106', 'Four waters', 'RMSD (Å)', 'Top pose ≤ 1.0 Å', 'Modes of the rotated pose', 'Score (kcal mol⁻¹)', 'Zn–N (Å)',
       'Strict gate', '2.0 Å criterion'],
      [':----------------------', ':----------', ':-----', '----------:', '--------:', '-------:', '-------:', '-----:', ':-----', ':-------'], rows,
      'Docking of sulfanilamide into carbonic anhydrase II (6RL9) with Glu106 as acid or carboxylate, without and with four crystal waters. '
      'The acid is the state of the receptors behind every docking number of the main text and of Section S13 for which the carboxylate is not named; the waters are 423, 516, 553 and 568 of 6RL9, kept as fixed receptor atoms with '
      'hydrogens placed for the largest number of hydrogen bonds. RMSD is the symmetry-corrected heavy-atom RMSD of the top pose from the crystal ligand (mean ± SD over the seeds), followed by the '
      'number of seeds whose top pose lies within 1.0 Å of it. A mode of the rotated pose lies beyond 2.0 Å from the crystal ligand and within 2.0 Å of the ring-rotated pose of the dry receptor '
      f"({v('dock.wat4.rot_ref')} of the first row); it is counted over every mode of every seed. The score is the mean score of the top poses, and Zn–N the mean distance from the zinc to their sulfonamide "
      'nitrogen. The gates of Section S13 apply to the validation runs (ETKDG conformers) only.', 's-dock-glu106-water',
      notes=['^a^ The hydrogen-bond ring of the four waters reversed. ^b^ The waters as bare oxygen atoms, as in the control receptor of Section S13. '
             f"Exhaustiveness {v('protocol.dock_exhaustiveness')}, {f(v('protocol.dock_box'), 0)} Å box, ligand inputs identical within each block of rows; AD4Zn, AutoDock4Zn. "
             f"The SD of the score over the seeds is at most {f(max(_sds['ad4zn']), 3)} kcal mol⁻¹ with AD4Zn, {f(max(_sds['vina']), 3)} with Vina and {f(max(_sds['vinardo']), 3)} with Vinardo. "
             'The rows with the waters are a diagnostic and not results of the protocol: a rigid receptor cannot let a ligand displace water, so the rotated pose is excluded by construction, '
             'the four waters were always kept together, and a pass of a gate in such a receptor does not validate the protocol. '
             f"Stored modes of the rotated pose scored without a search in the maps with the waters give +{_n(f(v('dock.wat4.fixed.rot_min_all'), 0))} to "
             f"+{_n(f(v('dock.wat4.fixed.rot_max_all'), 0))} kcal mol⁻¹ over every placement of the water hydrogens tried; in the receptors of the table the best near-native mode of each seed changes by "
             f"{f(v('dock.wat4.fixed.acid.near_change'), 3)} (acid) and {f(v('dock.wat4.fixed.carb.near_change'), 3)} kcal mol⁻¹ (carboxylate) on average."],
      doc='si', section='S15')

rows, _sds = [], []
for pdb_, lig_, kind_, what_ in (('3HS4', 'sfa_neutral', 'val', 'Sulfanilamide, ETKDG, cross-dock'), ('3HS4', 'sfa_neutral', 'prod', 'Sulfanilamide, DFT geometry, cross-dock'),
                                 ('6RL9', 'azm_neutral', 'val', 'Acetazolamide, ETKDG, cross-dock'), ('3HS4', 'azm_neutral', 'val', 'Acetazolamide, ETKDG, redock'),
                                 ('3HS4', 'azm_anion', 'val', 'Acetazolamide anion, ETKDG, redock')):
    for i_, glu_ in enumerate(('acid', 'carb')):
        k_ = f'dock.glu106.{pdb_}.{lig_}.ad4zn.{kind_}.{glu_}'
        _sds.append(v(k_ + '.score_sd'))
        if lig_ == 'sfa_neutral':            # every sulfanilamide top pose that is not near-native is the ring-rotated orientation
            assert v(k_ + '.top_between') == 0 and v(k_ + '.top_near') + v(k_ + '.top_beyond2') == v(k_ + '.n_seeds'), k_
        rows.append([pdb_ if i_ == 0 else '', f"{what_}, {v(k_ + '.n_seeds')} seeds" if i_ == 0 else '', _GLU[glu_], _pm(k_), f"{v(k_ + '.n_le_1')} of {v(k_ + '.n_seeds')}",
                     f(v(k_ + '.score'), 3, exact=True), f(v(k_ + '.zn_n'), 2), f"{v(k_ + '.zn_n_pass')} of {v(k_ + '.n_seeds')}"] + _gate2(k_))
_fn = [v(f'dock.glu106.fixed.{p_}.{q_}.near') for p_ in ('6RL9', '3HS4') for q_ in ('val', 'prod')]
_fr = [v(f'dock.glu106.fixed.{p_}.{q_}.rot') for p_ in ('6RL9', '3HS4') for q_ in ('val', 'prod')]
table(ST / 'tableS32_dock_glu106_other.md',
      ['Receptor', 'Ligand, input', 'Glu106', 'RMSD (Å)', 'Top pose ≤ 1.0 Å', 'Score (kcal mol⁻¹)', 'Zn–N (Å)', 'Seeds with Zn–N ≤ 2.3 Å',
       'Strict gate', '2.0 Å criterion'],
      [':-------', ':----------------------', ':----------', '----------:', '--------:', '-------:', '-----:', '--------:', ':-----', ':-------'], rows,
      'The other AutoDock4Zn runs of the carbonic anhydrase II receptors repeated with Glu106 as a carboxylate. Columns as in @tbl:s-dock-glu106-water. For sulfanilamide every top pose that '
      'is not within 1.0 Å of the reference lies beyond 2.0 Å, in the ring-rotated orientation. The reference of a cross-docking run is the crystal ligand of the other structure carried over by '
      'Cα superposition.', 's-dock-glu106-other',
      notes=['The acid rows are the runs of the main text and of Section S13; ligand inputs are identical within each pair of rows. '
             f"The SD of the score over the seeds is at most {f(max(_sds), 3)} kcal mol⁻¹. "
             f"Stored modes of the four sulfanilamide runs scored without a search in both sets of maps change by {f(max(_fn), 3)} to {f(min(_fn), 3)} kcal mol⁻¹ (near-native modes) and "
             f"{f(max(_fr), 3)} to {f(min(_fr), 3)} kcal mol⁻¹ (modes beyond 2.0 Å) with the carboxylate, so the ring-rotated orientation gains "
             f"{f(v('dock.glu106.fixed.rot_gain_min'), 3)} to {f(v('dock.glu106.fixed.rot_gain_max'), 3)} kcal mol⁻¹, the size of the separation of the two orientations. "
             f"The Vina and Vinardo scores of {v('dock.glu106.fixed.vina.n_poses')} stored top poses do not change ({f(v('dock.glu106.fixed.vina.max_change'), 3)} kcal mol⁻¹), and no verdict of the "
             f"{v('dock.glu106.n_runs_vina')} repeated Vina and Vinardo runs changes. Not repeated with the carboxylate: the six input-geometry runs of sulfanilamide, the 20-seed and "
             'higher-exhaustiveness variants of the 6RL9 validation run, every run of sulfamethoxazole, *p*-aminobenzoate and 4-sulfamoylanilinium, the other Vina validation variants, and the His64 variants. '
             'The heavy atoms of the receptor were not moved, and the test does not establish the protonation state of Glu106.'],
      doc='si', section='S15')

# S33 simulation settings
_kq = lambda i_: f"{f(v(f'md.prot.eq{i_}.length'), 0)} ps, {f(v(f'md.prot.eq{i_}.timestep'), 0)} fs"   # noqa: E731
_hk = lambda i_: (f"{v(f'md.prot.eq{i_}.held')}, {_n(v(f'md.prot.eq{i_}.k'))} kJ mol⁻¹ nm⁻²" if v(f'md.prot.eq{i_}.k') else 'nothing held')   # noqa: E731
_HELD = {'solute heavy': 'heavy atoms of the solute', 'protein backbone': 'protein backbone'}
rows = [
    ['*System*', ''],
    ['Structure', f"PDB {v('md.sys.pdb')}, chain A: {v('md.sys.residues')} residues (His{v('md.sys.res_first')} to Lys{v('md.sys.res_last')}, no residue {v('md.sys.res_missing')}), acetyl cap on His{v('md.sys.res_first')}; Glu106 as a carboxylate, His64 neutral (HIE)"],
    ['Ligand', f"sulfonamidate anion of sulfanilamide ({v('md.sys.ligand_atoms')} atoms), N2 bonded to the zinc; complex charge {f(v('md.sys.charge_solute'), 0)}"],
    ['Start X', f"crystal pose (d1 {f(v('md.sys.X.d1'), 2)} Å, d2 {f(v('md.sys.X.d2'), 2)} Å)"],
    ['Start N', f"docked near-native pose, {f(v('md.sys.N.dock_rmsd'), 2)} Å from the crystal ligand (d1 {f(v('md.sys.N.d1'), 2)} Å, d2 {f(v('md.sys.N.d2'), 2)} Å)"],
    ['Start R', f"docked ring-rotated pose, {f(v('md.sys.R.dock_rmsd'), 2)} Å from the crystal ligand (d1 {f(v('md.sys.R.d1'), 2)} Å, d2 {f(v('md.sys.R.d2'), 2)} Å)"],
    ['Atoms (X / N / R)', ' / '.join(_n(v(f'md.sys.{s_}.atoms')) for s_ in 'XNR')],
    ['Water molecules (X / N / R)', ' / '.join(_n(v(f'md.sys.{s_}.waters')) for s_ in 'XNR')],
    ['Crystal waters kept (X / N / R)', ' / '.join(str(v(f'md.sys.{s_}.xtal_waters')) for s_ in 'XNR') + f" of {v('md.sys.xtal_waters_in_entry')} sites"],
    ['Ions', f"{v('md.sys.na')} Na⁺ and {v('md.sys.cl')} Cl⁻ ({f(v('md.sys.salt_target'), 2)} M counted in the solvent volume)"],
    ['Cell', f"{v('md.sys.box')}, {f(v('md.sys.box_buffer'), 0)} Å buffer, edge {f(v('md.sys.box_edge'), 1)} Å as built"],
    ['*Force field*', ''],
    ['Protein, water, ions', f"{v('md.ff.protein')}, {v('md.ff.water')}, {v('md.ff.ions')}"],
    ['Ligand', f"{v('md.ff.ligand')}; equilibrium values of the bond S–N2 and of the angles C4–S–N2 and S–N2–H and the torsion about S–N2 (a series on H–N2–S–C4, H–N2–S–O set to zero) "
               'refitted to the DFT of this work'],
    ['Zinc site', f"bonded model (MCPB.py, AmberTools {v('md.sw.ambertools')}): {v('md.zn.2e.n_bonds')} Zn–N bonds, {v('md.zn.2e.n_angles')} angles, {v('md.zn.n_torsions')} torsions with zero force constant"],
    ['Zinc set 2e (nine runs)', f"equilibrium values of 6RL9 (Zn–N2 {f(v('md.zn.2e.bond.N2.r0'), 3)} Å; Zn–N(His) {f(v('md.zn.2e.bond.His94.r0'), 3)}, {f(v('md.zn.2e.bond.His96.r0'), 3)}, "
                                f"{f(v('md.zn.2e.bond.His119.r0'), 3)} Å); force constants {f(v('md.zn.2e.bond.N2.k'), 1)} and {f(v('md.zn.2e.bond.His94.k'), 1)}, {f(v('md.zn.2e.bond.His96.k'), 1)}, "
                                f"{f(v('md.zn.2e.bond.His119.k'), 1)} kcal mol⁻¹ Å⁻², {f(v('md.zn.2e.k_N_Zn_N'), 0)} (N–Zn–N) and {f(v('md.zn.2e.k_Zn_N_C'), 0)} kcal mol⁻¹ rad⁻² (Zn–N–C)"],
    ['Zinc set 2s (one run)', f"Seminario force constants and the equilibrium values of the optimised model (Zn–N2 {f(v('md.zn.2s.bond.N2.r0'), 4)} Å, {f(v('md.zn.2s.bond.N2.k'), 1)} kcal mol⁻¹ Å⁻²)"],
    ['Terms at N2, as built', f"angles Zn–N2–S ({f(v('md.zn.n2.Zn_N2_S.k'), 0)} kcal mol⁻¹ rad⁻², {f(v('md.zn.n2.Zn_N2_S.theta0'), 2)}°, the angle of 6RL9), Zn–N2–H ({f(v('md.zn.n2.Zn_N2_H.k'), 0)}, "
                              f"{f(v('md.zn.n2.Zn_N2_H.theta0'), 1)}°, the angle of the constructed hydrogen in start X) and S–N2–H ({f(v('md.ff.lig.angle_HN21_N2_S.k'), 2)}, {f(v('md.ff.lig.angle_HN21_N2_S.eq'), 1)}°), "
                              'with the torsion H–N2–S–C4 of the ligand; the same in both zinc sets'],
    ['Charges of the zinc site', f"RESP on the {v('md.sys.zn_large_atoms')}-atom capped model, {v('md.qm.method')}; Zn {f(v('md.zn.charge_zn'), 2, sign=True)} e"],
    ['*Dynamics*', ''],
    ['Program', f"OpenMM {v('md.sw.openmm')}, CUDA, {v('md.sw.precision')} precision"],
    ['Non-bonded', f"particle-mesh Ewald, cut-off {f(v('md.prot.cutoff'), 1)} nm, tolerance {f(v('md.prot.ewald_tol'), 4)}, dispersion correction"],
    ['Constraints, masses', f"bonds to hydrogen constrained, rigid water; hydrogen mass {f(v('md.prot.h_mass'), 3)} u"],
    ['Thermostat, barostat', f"Langevin (middle scheme), {f(v('md.prot.temperature'), 0)} K, {f(v('md.prot.friction'), 0)} ps⁻¹; Monte Carlo, {f(v('md.prot.pressure'), 0)} bar"],
    ['Minimisation', f"heavy atoms of the solute and crystal-water oxygens held, {_n(v('md.prot.min_k'))} kJ mol⁻¹ nm⁻²"],
    ['Stage 1, heating', f"{_kq(1)}, NVT, {f(v('md.prot.eq1.T_start'), 0)} to {f(v('md.prot.temperature'), 0)} K; {_HELD[v('md.prot.eq1.held')]} held, {_n(v('md.prot.eq1.k'))} kJ mol⁻¹ nm⁻²"],
    ['Stage 2', f"{_kq(2)}, NPT; {_HELD[v('md.prot.eq2.held')]} held, {_n(v('md.prot.eq2.k'))} kJ mol⁻¹ nm⁻²"],
    ['Stage 3', f"{_kq(3)}, NPT; {_HELD[v('md.prot.eq3.held')]} held, {_n(v('md.prot.eq3.k'))} kJ mol⁻¹ nm⁻²; extended by {f(v('md.prot.hyd.ext_length'), 0)} ps if no water lies in the place of the crystal ring (start R)"],
    ['Stage 4, ligand released', f"{_kq(4)}, NPT; {_HELD[v('md.prot.eq4.held')]} held, {_n(v('md.prot.eq4.k'))} kJ mol⁻¹ nm⁻²"],
    ['Stage 5', f"{_kq(5)}, NPT; {_HELD[v('md.prot.eq5.held')]} held, {_n(v('md.prot.eq5.k'))} kJ mol⁻¹ nm⁻²"],
    ['Stage 6', f"{_kq(6)}, NPT; nothing held"],
    ['Production', f"{f(v('md.prot.prod.length'), 0)} ns, NPT, {f(v('md.prot.prod.timestep'), 0)} fs, no restraint and no bias"],
    ['Frames written', f"solute every {f(v('md.prot.out.solute'), 0)} ps ({_n(v('md.prot.frames_solute'))} per run), all atoms every {f(v('md.prot.out.full'), 0)} ps ({_n(v('md.prot.frames_full'))} per run)"],
    ['Runs', f"{v('md.prot.n_runs_primary')} with set 2e ({v('md.prot.n_runs_per_start')} per start), {v('md.prot.n_runs_sensitivity')} sensitivity runs from start R, {v('md.prot.n_runs_validation')} validation runs of "
             f"{f(v('md.prot.validation_ns'), 0)} ns; {_n(v('md.prot.total_production_ns'))} ns of production in all"],
    ['*Analysis*', ''],
    ['Orientation', f"d1, N1 to Phe131 CZ; d2, N1 to Asn67 OD1. Crystal-like: d1 < {f(v('md.def.mid_d1'), 1)} Å and d2 > {f(v('md.def.mid_d2'), 1)} Å; rotated-like: d1 > {f(v('md.def.mid_d1'), 1)} Å and "
                    f"d2 < {f(v('md.def.mid_d2'), 1)} Å; intermediate otherwise"],
    ['Verdict of a run', f"CONVERTED or STAYED if at least {f(v('md.def.verdict_pct'), 0)} % of the frames of the last {f(v('md.def.window'), 0)} ns are in the state opposite to the start or in that of the start; MIXED otherwise"],
    ['Time of a conversion', f"first frame of the first visit of at least {f(v('md.def.dwell'), 0)} ns in the opposite state, counted from the release of the ligand"],
    ['Ligand RMSD', f"{v('md.sys.ligand_heavy_atoms')} heavy atoms against the crystal pose after the fit of the {v('md.sys.residues')} Cα atoms; smallest of the four equivalent labellings"],
    ['Hydrogen bond to Thr199', f"N2···OG1 at most {f(v('md.def.hb_dist'), 1)} Å and N2–H···OG1 at least {f(v('md.def.hb_angle'), 0)}°"],
    ['Place of a crystal water', f"occupied in an all-atom frame when a water oxygen lies within {f(v('md.def.water_site'), 1)} Å of the crystal position after the fit of the Cα atoms"],
    ['His64 deep in the site', f"Zn–NE2(His64) at most {f(v('md.def.his64_deep'), 1)} Å with χ1 between 0° and {f(v('md.def.his64_chi1_max'), 0)}°; the limit was chosen after looking at the data"],
    ['Nearly planar nitrogen', f"signed volume of the unit vectors from N2 to Zn, H and S below {f(v('md.def.n2_flat'), 1)} in size"],
]
table(ST / 'tableS33_md_settings.md', ['Quantity', 'Value'], [':------------------------------', ':-------------------------------------------------------------------------------------------'], rows,
      'Settings of the molecular dynamics of the carbonic anhydrase II–sulfanilamide complex. N1 is the aniline nitrogen and N2 the sulfonamide nitrogen of the ligand. '
      'Restraints are harmonic position restraints on the coordinates of the built system (*U* = *k* *d*²/2). '
      'The force constants of bonds and angles are those of the Amber energy function, *U* = *k*(*x* − *x*~0~)², which has no factor ½.', 's-md-settings',
      notes=['The equilibrium Zn–N distances and N–Zn–N angles of set 2e are those of 6RL9, so the zinc geometry of the nine runs is put in and is no result. '
             'The force-field terms at N2 failed two quantum-chemical tests (@tbl:s-md-gates) and are used as built. '
             f"The cell holds {f(v('md.sys.salt_in_box'), 2)} M NaCl when the salt is counted in its whole volume."],
      doc='si', section='S16')

# S34 the eleven runs
rows = []
for rid_ in ('X1', 'X2', 'X3', 'N1', 'N2', 'N3', 'R1', 'R2', 'R3', 'R1_2s', 'R1_n2fit'):
    b_ = f'md.run.{rid_}'
    rule_ = (f"{v(b_ + '.rule.ps')}; {v(b_ + '.rule.stage')}" + ('^c^' if v(b_ + '.verdict') == 'MIXED' else '')) if b_ + '.rule.ps' in F else 'none'
    sign_ = f"{v(b_ + '.sign.ps')}; {v(b_ + '.left.ps')}" if b_ + '.sign.ps' in F else '–'
    if b_ + '.sign.ps' not in F:
        assert v(b_ + '.n2.n_changes') == 0 and v(b_ + '.whole.n_R') == 0, rid_
    rows.append([{'R1_2s': 'R1^a^', 'R1_n2fit': 'R1^b^'}.get(rid_, rid_), v(b_ + '.start'), v(b_ + '.zinc_set'), v(b_ + '.n2_terms').split(',')[0], v(b_ + '.verdict'),
                 ' / '.join(f(v(f'{b_}.last20.pct_{s_}'), 2) for s_ in 'CRI'), rule_, sign_, f"{f(v(b_ + '.rmsd.mean'), 2)} ± {f(v(b_ + '.rmsd.sd'), 2)}", f(v(b_ + '.rmsd.last_half'), 2),
                 f(v(b_ + '.hb.pct'), 2), f(v(b_ + '.ca_rmsd.mean'), 2)])
table(ST / 'tableS34_md_runs.md',
      ['Run', 'Start', 'Zinc set', 'Terms at N2', 'Verdict', 'Last 20 ns: C / R / I (%)', 'Time of the rule (ps); stage', 'N2 inverts; ring leaves (ps)', 'Ligand RMSD (Å)',
       'Last half (Å)', 'H-bond to Thr199 (%)', 'Cα RMSD (Å)'],
      [':---', ':---', ':---', ':-------', ':------------', '---------------:', '------:', '-------:', '----------:', '----:', '------:', '-----:'], rows,
      'The eleven molecular dynamics runs of 100 ns. Start X is the crystal pose, N the docked near-native pose and R the docked ring-rotated pose. '
      # 11 Oct 2026 (repair round, E8/C11): the run names N1 and N2 are not the atom labels.
      'The rows N1 to N3 are the runs from start N; in the column headings N2 is the sulfonamide nitrogen of the ligand. '
      'C, R and I are the shares of crystal-like, rotated-like and intermediate frames of the last 20 ns, which give the verdict. '
      'The time of the rule is the first frame of the first visit of at least 1 ns in the crystal-like state, counted from the release of the ligand, with the equilibration stage in which it lies; '
      'the next column gives, for the runs from start R, the first frame with the configuration of the crystal pose at N2 and the first frame that is not rotated-like '
      '(the runs from X and N keep the configuration of the crystal pose and have no rotated-like frame). '
      'Ligand RMSD: heavy atoms against the crystal pose of 6RL9, mean ± SD over the production and mean of its last half. Hydrogen bond N2–H···OG1(Thr199): share of the production frames. '
      'Cα RMSD from 6RL9: mean over the production.', 's-md-runs',
      notes=['^a^ Sensitivity run with zinc set 2s, which failed its validation run in one part (@tbl:s-md-gates). ^b^ Sensitivity run with candidate terms at N2 that were not accepted '
             '(@tbl:s-md-gates). ^c^ The verdict is MIXED, so this is no conversion but the start of the first crystal-like visit of at least 1 ns. '
             f"All times of the table lie in the equilibration, before the production: the ligand is released at the start of stage 4 and the production starts {_n(v('md.prot.eq_released'))} ps later. "
             f"In stage 4 the protein backbone is held with {_n(v('md.prot.eq4.k'))} kJ mol⁻¹ nm⁻², in stage 5 with {_n(v('md.prot.eq5.k'))}, in stage 6 nothing is held; in all five runs from start R the "
             f"first frame with the configuration of the crystal pose at N2 and the first frame that is not rotated-like lie in stage {v('md.grp.R.event_stage')}. "
             f"The production of the nine runs with set 2e holds {v('md.grp.nine.n_R_frames')} rotated-like frames of {_n(v('md.grp.nine.n_frames'))}. "
             'The shares are not populations, and the times are single events with frames 10 ps apart: no rate, barrier or free energy follows from them. '
             'That the crystal-like state is kept and not left again belongs to the terms at N2 as built, which hold the nitrogen too stiffly; the two sensitivity runs enter no count.'],
      doc='si', section='S16')

# S35 the tests of the force field
_d = lambda k_, n_=2, **kw: f(v(k_), n_, **kw)   # noqa: E731
_gz = lambda t_: (f"Zn–N2 {_d(f'md.gate.GZ.{t_}.zn_n2', 3)} Å; Zn–N(His) {_d(f'md.gate.GZ.{t_}.zn_h94', 3)}, {_d(f'md.gate.GZ.{t_}.zn_h96', 3)}, {_d(f'md.gate.GZ.{t_}.zn_h119', 3)} Å; "   # noqa: E731
                  f"largest angle departure {_d(f'md.gate.GZ.{t_}.ang_max_dev', 1, sign=True)}°; oxygen within 2.6 Å in {v(f'md.gate.GZ.{t_}.so_n')} of {v('md.gate.GZ.n_frames')} frames")
_gvv = lambda t_: (f"{_d(f'md.gate.GV.{t_}.rmsd')} Å; {_d(f'md.gate.GV.{t_}.hb_pct', 1)} %; Zn–N2 {_d(f'md.gate.GV.{t_}.zn_n2', 3)} Å, largest angle departure "   # noqa: E731
                   f"{_d(f'md.gate.GV.{t_}.ang_max_dev', 1, sign=True)}°; oxygen within 2.6 Å in {v(f'md.gate.GV.{t_}.any_n')} of {v('md.gate.GV.n_full')} all-atom frames ({_d(f'md.gate.GV.{t_}.any_pct', 0)} %); "
                   f"{_d(f'md.gate.GV.{t_}.c4sn2', 1)}° and {_d(f'md.gate.GV.{t_}.tau_fold', 1)}°")
_win = lambda k_: f"{f(v(k_)[0], 2)} to {f(v(k_)[1], 2)}"   # noqa: E731
_GZC = (f"Zn–N2 {_win('md.gate.GZ.limit.zn_n2')} Å; Zn–N(His) {_win('md.gate.GZ.limit.zn_his')} Å; each mean N–Zn–N angle within "
        f"{_d('md.gate.GZ.limit.angle', 0)}° of 6RL9; an oxygen within 2.6 Å of the zinc in at most {_d('md.gate.GZ.limit.fifth_pct', 0)} % of the frames")
rows = [
    ['G-L', 'Ligand terms at sulfur against the DFT of this work (sulfonamidate, implicit water)',
     f"C–S–N within {_d('md.gate.GL.limit.csn', 0)}° and S–N within {_d('md.gate.GL.limit.sn', 2)} Å at the minimum; torsion profile within {_d('md.gate.GL.limit.scan', 1)} kcal mol⁻¹",
     f"{_d('md.gate.GL.after.csn', 2, sign=True)}°, {_d('md.gate.GL.after.sn', 3, sign=True)} Å, {_d('md.gate.GL.after.scan')} kcal mol⁻¹ (GAFF2 as assigned: {_d('md.gate.GL.before.csn', 2, sign=True)}°, "
     f"{_d('md.gate.GL.before.sn', 3, sign=True)} Å, {_d('md.gate.GL.before.scan')})", v('md.gate.GL.after.outcome') + ' after the refit'],
    ['G-Z, set 2e', 'Zinc site in the first 100 ps of the validation run, with the site still held', _GZC, _gz('2e'), v('md.gate.GZ.2e.outcome')],
    ['G-Z, set 2s', 'as above', 'as above', _gz('2s'), v('md.gate.GZ.2s.outcome')],
    ['G-V, set 2e', f"Validation run of {f(v('md.prot.validation_ns'), 0)} ns from the crystal pose",
     f"ligand RMSD at most {_d('md.gate.GV.limit.rmsd', 1)} Å; hydrogen bond to Thr199 in at least {_d('md.gate.GV.limit.hb_pct', 0)} % of the frames; zinc site as in G-Z; C4–S–N2 within "
     f"{_d('md.gate.GV.limit.c4sn2', 0)}° of {_d('md.xtal.c4sn2', 1)}° and *τ*~fold~ within {_d('md.gate.GV.limit.tau_fold', 0)}° of {_d('md.xtal.tau_fold', 1)}°", _gvv('2e'), v('md.gate.GV.2e.outcome')],
    ['G-V, set 2s', 'as above', 'as above', _gvv('2s'), v('md.gate.GV.2s.outcome') + ' (the oxygen count)'],
    ['G-R', f"Terms at N2 against {v('md.qm.gate_method')} energies of the {v('md.sys.zn_small_atoms')}-atom zinc-site model, {v('md.gate.GR.n')} comparisons",
     f"root mean square of force field minus quantum chemistry at most {_d('md.gate.GR.limit', 1)} kcal mol⁻¹",
     f"{_d('md.gate.GR.rms')} kcal mol⁻¹ ({_d('md.gate.GR.rms_no_d3')} without D3). Force field against quantum chemistry, relative to the crystal-like placement: planar nitrogen "
     f"{_d('md.gate.GR.planar_N.mm')} against {_d('md.gate.GR.planar_N.qm')}; hydrogen moved into the plane {_d('md.gate.GR.H_in_plane.mm')} against {_d('md.gate.GR.H_in_plane.qm')}; ligand turned by 30° about "
     f"Zn–N2 {_d('md.gate.GR.turn_30.mm')} against {_d('md.gate.GR.turn_30.qm')}; rotated placement {_d('md.gate.GR.Rx.mm')} against {_d('md.gate.GR.Rx.qm')} ({_d('md.gate.GR.Rx.mm_minus_qm', 2, sign=True)})",
     v('md.gate.GR.outcome')],
    ['G-R, refit', f"Refit of the terms at N2: {v('md.gate.GR2.n_candidates')} candidates fitted to {v('md.gate.GR2.n_fit_rows')} single-point energies ({v('md.gate.GR2.n_fit_old')} of G-R and "
                   f"{v('md.gate.GR2.n_fit')} new ones)",
     f"leave-one-out error at most {_d('md.gate.GR2.limit.loo', 1)} kcal mol⁻¹; angles at N2 of the force-field minimum within {_d('md.gate.GR2.limit.angle', 0)}° of the quantum-chemical optimum",
     f"best leave-one-out errors {_d('md.gate.GR2.loo.as_built_charges')} and {_d('md.gate.GR2.loo.two_placements')} kcal mol⁻¹ (charges as built; charges fitted to both placements); "
     f"Zn–N2–S of the force-field minimum off by {_d('md.gate.GR2.angle_miss', 1)}° and {_d('md.gate.GR2.angle_miss_two_placements', 1)}°",
     v('md.gate.GR2.outcome')],
]
table(ST / 'tableS35_md_gates.md', ['Test', 'What is tested', 'Criterion', 'Found', 'Outcome'],
      [':--------', ':----------------------', ':---------------------------------', ':-----------------------------------------------', ':-------------'], rows,
      'Tests of the force field of the simulations, with those that failed. Energies in kcal mol⁻¹ are relative to the crystal-like placement of the ligand in the bare zinc-site model '
      '(zinc, three ethylimidazoles and the ligand, in vacuum). Zinc set 2e is the set of the nine runs and set 2s that of one sensitivity run.', 's-md-gates',
      notes=[f"G-R failed, and no candidate of its refit was accepted, so the nine runs use the terms at N2 as built; the {v('md.gate.GR2.n_test')} points of the test set of the refit were never opened, "
             'because the verdict was fixed without them. Direction of the errors of G-R: the force field resists a planar nitrogen about twice as strongly as quantum chemistry, makes the turn about '
             f"Zn–N2 too easy, and puts the rotated placement {_d('md.gate.GR.Rx.mm_minus_qm')} kcal mol⁻¹ too high. The last error has the sign of the outcome of the runs from start R; its sign rests on the "
             f"dispersion term of the reference ({_d('md.gate.GR.Rx.mm_minus_qm_no_d3', 2, sign=True)} without D3), it falls to {_d('md.gate.GR.Rx.mismatch_two_placements', 2, sign=True)} with charges fitted "
             'to both placements, and its size in the protein is not known. '
             f"Set 2s missed G-V by one all-atom frame ({v('md.gate.GV.2s.any_n')} against the {v('md.gate.GV.n_full_allowed')} allowed); it is not a five-coordinate model. "
             'For set 2e the zinc parts of G-Z and G-V show only that the rest of the force field does not move the site away from the values put in. '
             'The zinc-site model of G-R and of its refit has no Thr199, no pocket wall and no water. The files of the supplementary archive name the tests as here and the refit `G-R2`.'],
      doc='si', section='S16')

# ------------------------------------------------------------------------------------------------ index
L = ['# Table index (generated by paper/main/data/build_tables.py)', '',
     'Include a table in a manuscript with `<!-- include: tables/<file> -->` where the table belongs; the label '
     'is set inside each file. Numbers come from `paper/main/data/facts.json` and the files it cites.', '',
     '| Document | File | Label | Section | Caption (first sentence) |', '|---|---|---|---|---|']
for doc, fn, lab, sec, cap in INDEX:
    L.append(f'| {doc} | `{fn}` | `tbl:{lab}` | {sec} | {cap} |')
L += ['', 'SI notes: section S2 also has `coordinates_S2.md` (coordinate blocks, not a table); S8 (ELF/LOL) is figures only; '
      'S14 (ADMET) is Table S_admet and the condensed Fukui functions of S9 are Table S_fukui (added 29 Sep 2026). The SI skeleton has no dedicated section for charges, MEP, QTAIM/IGMH and interaction '
      'energies; they are placed under S9 here. Printed SI numbers follow the order of the include lines in `paper/si/manuscript.md`, '
      'not the file names: after VM burst 3 (29 Sep 2026) `tableS13b_periodic_profile.md` prints as Table S15, `tableS18b_cluster_qtaim.md` '
      'as Table S22 and `tableS26b_torsion_scan.md` as Table S30, and the later files shift accordingly (see the header of the SI manuscript).']
(MT / 'README.md').write_text('\n'.join(L) + '\n')
(ST / 'README.md').write_text('\n'.join(L) + '\n')
print(f'{sum(1 for x in INDEX if x[0] == "main")} main tables, {sum(1 for x in INDEX if x[0] == "si")} SI tables')
