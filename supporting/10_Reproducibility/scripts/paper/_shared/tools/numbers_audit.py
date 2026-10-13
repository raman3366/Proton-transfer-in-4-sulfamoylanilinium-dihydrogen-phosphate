#!/usr/bin/env python3
"""numbers_audit.py -- every number in a manuscript must come from facts.json or from the literature notes.

Usage
  python3 paper/_shared/tools/numbers_audit.py <main | si | dir> [--manuscript FILE] [--facts FILE] [--lit DIR]
          [--no-auto-diff] [--quiet] [--json out.json]

  main / si  ->  paper/<main|si>/manuscript.md, checked against paper/main/data/facts.json (one fact file for both
  documents) and the literature notes in lit/ (REVIEW.md, notes_*.md). A directory holding manuscript.md works
  too (facts and notes default to the same shared files). Adapted from cqc-dft (one directory per paper).

What is extracted
  Every number in the title, abstract, keywords, body text, figure and table captions, table cells and back
  matter, with its unit and 40 characters of context on each side. Not extracted: the reference list,
  display equations (and integers inside inline math), figure/table/equation/section/reference numbers
  (Fig. 2a, Figs. 1–3, Table S1, Eq. (4), Section 3.1, heading numbers), years inside citations (a year that is
  part of a cited key, or after "et al."), atom labels and formulae (N4, O11, C20H15N7O7, H₂O), basis-set and
  functional names (6-311++G(d,p), M06-2X, B3LYP-D3), locants (1,2,4-triazole, 3-amino, benzimidazol-3-ium),
  standard uncertainties "(3)", ratios (1:1:1), subscripts/superscripts written ~2~ or ^2^, DOIs.

How a number is matched (first rule that applies)
  1. facts.json: key -> {"value", "unit", "case", "source", "note"}. Numeric values, the numeric members of
     dict/list values (key.member, key[i]) and numbers inside string values (e.g. "CCDC 2265990") are
     candidates. A printed number x with d decimals matches a value v when |x| and |v| agree to within half a
     unit of the last printed digit (rounding to the printed precision; sign or absolute value are both
     accepted, and a sign flip is shown). After "about", "≈", "~", "nearly", ... an integer ending in zeros is
     compared at its significant figures (about 3500 -> ±50).
     Units: the printed unit (after the number, inherited across "1.05 and 1.70 Å", "−0.07 to +0.09 a.u.",
     "0.03 ± 0.01 Å", or taken from a table's column header / caption) must be compatible with the fact's:
     Ha ↔ eV (27.211386) ↔ kcal mol⁻¹ (627.5095 per Ha; 23.0605 per eV) ↔ kJ mol⁻¹, mHa, µHa; a.u. of energy
     or potential = Ha (a printed "a.u." is never a converted eV value); Å ↔ pm ↔ nm ↔ Bohr. Printed without a
     unit -> compared with every fact in its own unit. Printed with a unit -> never matched to a fact that has no
     unit (numbers inside string facts excepted). Table units come from the column header ("D–H (Å)"), from a
     caption exception naming the column symbol ("(eV; *S* in eV⁻¹)") or from the caption's single unit.
     When several facts fit, the one with the right unit, model (gas phase / IEF-PCM / X-ray words in the context),
     key words in the context and the closest value is shown, with the number of other candidates; a one-digit
     number with several candidates is marked "weak" so that a reviewer checks it.
  2. derived values listed in facts.json, either as a top-level "_derived" entry
         "_derived": {"pcm_minus_gas_gap": "pcm.gap_ev - gas.gap_ev", ...}      (or a list of expressions,
                                                                                 or of {"key","expr","unit"})
     or as an ordinary entry with "expr": "a - b" (or "derived_from": ["a", "b"] meaning a - b) and "unit".
     Expressions may use + - * / ( ) abs() numbers and fact keys.
     Automatic counterpart differences are added too (disable with --no-auto-diff): for facts whose keys differ
     only by the model prefix (gas., pcm., xray., periodic., saddle.) and that share a unit, prefix1.x - prefix2.x.
  3. literature: the number sits in a sentence (or table row/caption) that cites a literature key, and the same
     number (to the printed precision) appears in lit/REVIEW.md or lit/notes_*.md, preferably in a
     paragraph that names the cited key ("beside key"), otherwise anywhere in those files ("elsewhere in notes").
  4. otherwise UNMATCHED -> exit status 1.
"""
import argparse
import ast
import bisect
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mdtext import Doc, PAPERS, PAPER_IDS, LIT, MATH, rel  # noqa: E402

CASE_PREFIXES = ('gas', 'pcm', 'xray', 'periodic', 'saddle')

# ---------------------------------------------------------------------------------------------------------------
# units
# ---------------------------------------------------------------------------------------------------------------
HA_EV, HA_KCAL, HA_KJ, BOHR_A = 27.211386, 627.5095, 2625.4996, 0.529177210903
GROUPS = {
    'energy': {'Ha': 1.0, 'a.u.': 1.0, 'eV': 1 / HA_EV, 'kcal/mol': 1 / HA_KCAL, 'kJ/mol': 1 / HA_KJ,
               'mHa': 1e-3, 'µHa': 1e-6},
    'length': {'Å': 1.0, 'pm': 0.01, 'nm': 10.0, 'bohr': BOHR_A},
    'mass': {'g/mol': 1.0, 'Da': 1.0, 'amu': 1.0},
    'time': {'s': 1.0, 'min': 60.0, 'h': 3600.0},
}
UNIT_GROUP = {u: (g, f) for g, d in GROUPS.items() for u, f in d.items()}

FACT_UNIT_MAP = {
    'a': 'Å', 'å': 'Å', 'angstrom': 'Å', 'pm': 'pm', 'nm': 'nm', 'bohr': 'bohr',
    'ev': 'eV', 'ha': 'Ha', 'hartree': 'Ha', 'eh': 'Ha', 'mha': 'mHa', 'microha': 'µHa', 'µha': 'µHa', 'uha': 'µHa',
    'kcal mol-1': 'kcal/mol', 'kcal/mol': 'kcal/mol', 'kj mol-1': 'kJ/mol', 'kj/mol': 'kJ/mol',
    'cm-1': 'cm-1', 'deg': '°', '°': '°', 'degree': '°', 'degrees': '°', 'd': 'D', 'debye': 'D', 'e': 'e',
    'a.u.': 'a.u.', 'au': 'a.u.', '%': '%', 'k': 'K', 'ev-1': 'eV-1', 'a2': 'Å2', 'å2': 'Å2', 'å^2': 'Å2',
    'a3': 'Å3', 'å3': 'Å3', 'å^3': 'Å3', 'e a-3': 'e/Å3', 'e å-3': 'e/Å3', 'g cm-3': 'g/cm3', 'g mol-1': 'g/mol',
    'da': 'Da', 'amu': 'amu', 'min': 'min', 'h': 'h', 's': 's', 'pa': 'Pa', 'mm': 'mm', 'mm-1': 'mm-1',
    'ha/bohr': 'Ha/bohr', 'ha bohr-1': 'Ha/bohr', 'ha/bohr2': 'Ha/bohr2', 'ha bohr-2': 'Ha/bohr2',
    'kj mol-1 a-1': 'kJ/mol/Å',
}

# printed units, tried in this order right after a number
PRINTED_UNITS = [
    (r'kcal[\s  ]*(?:mol⁻¹|mol\^?-1\^?|mol−1|/[\s]*mol)', 'kcal/mol'),
    # 29 Sep 2026 (SAPH, VM burst 3): energy slopes of the periodic profile, kJ mol⁻¹ Å⁻¹ (facts unit 'kJ mol-1 A-1');
    # tried before kJ/mol so that the Å⁻¹ is not dropped
    (r'kJ[\s  ]*(?:mol⁻¹|mol\^?-1\^?|mol−1)[\s  ]*Å(?:⁻¹|\^?-1\^?|−1)', 'kJ/mol/Å'),
    (r'kJ[\s  ]*(?:mol⁻¹|mol\^?-1\^?|mol−1|/[\s]*mol)', 'kJ/mol'),
    (r'cm(?:⁻¹|\^?-1\^?|−1)', 'cm-1'),
    (r'eV(?:⁻¹|\^?-1\^?|−1)', 'eV-1'),
    (r'e[\s  ]*Å(?:⁻³|\^?-3\^?|−3)', 'e/Å3'),
    (r'g[\s  ]*cm(?:⁻³|\^?-3\^?|−3)', 'g/cm3'),
    (r'g[\s  ]*mol(?:⁻¹|\^?-1\^?|−1)', 'g/mol'),
    (r'mm(?:⁻¹|\^?-1\^?|−1)', 'mm-1'),
    (r'eV\b', 'eV'),
    (r'Å(?:²|\^2\^|2\b)', 'Å2'), (r'Å(?:³|\^3\^|3\b)', 'Å3'), (r'Å', 'Å'), (r'[Aa]ngstr[oö]ms?\b', 'Å'),
    (r'A(?![\w\'′])', 'Å'),
    (r'pm\b', 'pm'), (r'nm\b', 'nm'), (r'[Bb]ohr\b|a₀|a~0~', 'bohr'),
    (r'Ha[\s  ]*Bohr(?:⁻¹|\^?-1\^?|−1)|Ha/Bohr\b', 'Ha/bohr'),          # force (gradient) units
    (r'Ha[\s  ]*Bohr(?:⁻²|\^?-2\^?|−2)|Ha/Bohr(?:²|2)', 'Ha/bohr2'),   # force-constant (Hessian) units
    (r'mHa\b|mE~h~', 'mHa'), (r'[µμ]Ha\b', 'µHa'), (r'Ha\b|[Hh]artrees?\b|E~h~|Eh\b', 'Ha'),
    (r'a\.u\.|au\b', 'a.u.'),
    (r'Debye\b|D(?![\w\'′])', 'D'),
    (r'°C', '°C'), (r'°|deg\b|degrees?\b', '°'),
    (r'K(?![\w\'′])', 'K'), (r'%', '%'),
    (r'e(?![\w\'′])', 'e'),
    (r'Da\b', 'Da'), (r'amu\b', 'amu'), (r'min\b', 'min'), (r'h\b', 'h'), (r's\b', 's'),
    (r'Pa\b', 'Pa'), (r'GPa\b', 'GPa'), (r'MPa\b', 'MPa'), (r'mm\b', 'mm'), (r'ppm\b', 'ppm'),
]
PRINTED_UNITS = [(re.compile(r'[\s   ]?(?:' + rx + ')'), u) for rx, u in PRINTED_UNITS]


def fact_unit(u, member=None):
    if u is None:
        return None
    s = str(u).strip()
    if not s or s.lower() in ('dimensionless', 'none', 'count', 'points'):
        return None
    if ',' in s:                                       # composite unit of a dict value ("Å, deg")
        parts = [fact_unit(p) for p in s.split(',')]
        if member is not None and re.search(r'angle|deg|theta|phi|tors|dihed', member, re.I) and '°' in parts:
            return '°'
        if member is not None and re.search(r'D-H|H\.\.\.A|D\.\.\.A|dist|length|^r_|^d_', member):
            return 'Å' if 'Å' in parts else None
        return None
    s2 = re.sub(r'\s+', ' ', s.replace('^', '').replace('−', '-')).lower()
    s3 = re.sub(r'\s*\([^()]*\)$', '', s2)            # "a.u. (RDG)" -> "a.u."
    return FACT_UNIT_MAP.get(s2) or FACT_UNIT_MAP.get(s3) or s


def convert(v, uf, up):
    """value v in fact unit uf -> printed unit up; None if incompatible. (value, label)"""
    if uf == up:
        return v, ''
    if up == 'a.u.' and uf != 'Ha':
        return None, None                  # a printed "a.u." is a density or potential, not a converted eV value
    gf, gp = UNIT_GROUP.get(uf), UNIT_GROUP.get(up)
    if gf and gp and gf[0] == gp[0]:
        return v * gf[1] / gp[1], f'{uf}→{up}'
    return None, None


# ---------------------------------------------------------------------------------------------------------------
# facts
# ---------------------------------------------------------------------------------------------------------------
class Cand:
    __slots__ = ('key', 'value', 'unit', 'case', 'kind', 'decimals', 'text')

    def __init__(self, key, value, unit, case, kind, decimals=None, text=None):
        self.key, self.value, self.unit, self.case, self.kind = key, value, unit, case, kind
        self.decimals, self.text = decimals, text


def case_of(key, case):
    c = (str(case) if case else '').lower()
    # SAPH case names carry the environment as a token: sfa_pair_gas_d3, sfa_pair_pcm_b3lyp, cation_pcm_d3
    tokens = set(re.split(r'[_.\-]', c))
    for p in CASE_PREFIXES:
        if key.startswith(p + '.') or c == p or p in tokens:
            return p
    if 'saddle' in c:
        return 'saddle'
    return None


NUM_IN_STR = re.compile(r'(?<![\w.])[−\-]?\d+(?:\.\d+)?(?![\w.]*\d)')


class SafeExpr:
    OPS = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Add, ast.Sub, ast.Mult, ast.Div, ast.USub, ast.UAdd,
           ast.Constant, ast.Name, ast.Load, ast.Call)

    @classmethod
    def evaluate(cls, expr, values):
        keys = sorted(values, key=len, reverse=True)
        names = {}

        def sub(m):
            k = m.group()
            if k in values:
                n = f'_k{len(names)}'
                names[n] = values[k]
                return n
            return k
        e = re.sub(r'[A-Za-z_][\w.\[\]]*', sub, expr)
        tree = ast.parse(e, mode='eval')
        for node in ast.walk(tree):
            if not isinstance(node, cls.OPS):
                raise ValueError(f'unsupported syntax in {expr!r}')
            if isinstance(node, ast.Call) and not (isinstance(node.func, ast.Name) and node.func.id == 'abs'):
                raise ValueError(f'only abs() is allowed in {expr!r}')
            if isinstance(node, ast.Name) and node.id not in names and node.id != 'abs':
                raise ValueError(f'unknown key {node.id!r} in {expr!r}')
        return float(eval(compile(tree, '<expr>', 'eval'), {'__builtins__': {}, 'abs': abs}, names))


def load_facts(path, auto_diff=True):
    raw = json.loads(Path(path).read_text(encoding='utf-8'))
    cands, values, units, cases = [], {}, {}, {}
    errors = []
    derived_specs = []

    def add_num(key, v, unit, case, kind):
        if isinstance(v, bool) or not isinstance(v, (int, float)) or (isinstance(v, float) and not math.isfinite(v)):
            return
        cands.append(Cand(key, float(v), unit, case, kind))
        values[key] = float(v)
        units[key] = unit
        cases[key] = case

    for key, ent in raw.items():
        if key.startswith('_'):
            if key == '_derived':
                if isinstance(ent, dict):
                    derived_specs += [{'key': k, 'expr': v} if isinstance(v, str) else dict(v, key=k)
                                      for k, v in ent.items()]
                elif isinstance(ent, list):
                    derived_specs += [{'key': f'_derived[{i}]', 'expr': v} if isinstance(v, str) else v
                                      for i, v in enumerate(ent)]
            continue
        if not isinstance(ent, dict):
            ent = {'value': ent}
        v = ent.get('value')
        c = case_of(key, ent.get('case'))
        if ent.get('expr') or ent.get('derived_from'):
            expr = ent.get('expr') or ' - '.join(ent['derived_from'])
            derived_specs.append({'key': key, 'expr': expr, 'unit': ent.get('unit'), 'value': v})
            if v is None:
                continue
        if isinstance(v, dict):
            for mk, mv in v.items():
                add_num(f'{key}.{mk}', mv, fact_unit(ent.get('unit'), mk), c, 'member')
        elif isinstance(v, list):
            for i, mv in enumerate(v):
                add_num(f'{key}[{i}]', mv, fact_unit(ent.get('unit')), c, 'member')
        elif isinstance(v, str):
            for m in NUM_IN_STR.finditer(v):
                s = m.group().replace('−', '-')
                d = len(s.split('.')[1]) if '.' in s else 0
                cands.append(Cand(key, abs(float(s)), fact_unit(ent.get('unit')), c, 'string', d, v))
            for m in re.finditer(r'\d+(?:\.\d+){2,}', v):          # version strings
                cands.append(Cand(key, None, None, c, 'version', None, m.group()))
        else:
            add_num(key, v, fact_unit(ent.get('unit')), c, 'value')
    n_listed = 0
    for spec in derived_specs:
        try:
            val = SafeExpr.evaluate(spec['expr'], values)
        except (ValueError, SyntaxError, ZeroDivisionError) as e:
            errors.append(f'derived {spec.get("key")}: {e}')
            continue
        unit = fact_unit(spec.get('unit')) if spec.get('unit') is not None else None
        if unit is None:
            names = [k for k in values if re.search(r'(?<![\w.])' + re.escape(k) + r'(?![\w.])', spec['expr'])]
            us = {units[k] for k in names}
            unit = us.pop() if len(us) == 1 else None
        cands.append(Cand(f'{spec.get("key")} = {spec["expr"]}', val, unit, None, 'derived'))
        n_listed += 1
    n_auto = 0
    if auto_diff:
        groups = defaultdict(dict)
        for k in values:
            p, _, rest = k.partition('.')
            if p in CASE_PREFIXES and rest:
                groups[rest][p] = k
        for rest, byp in groups.items():
            ps = sorted(byp)
            for i in range(len(ps)):
                for j in range(i + 1, len(ps)):
                    ka, kb = byp[ps[i]], byp[ps[j]]
                    if units[ka] != units[kb]:
                        continue
                    cands.append(Cand(f'{ka} − {kb}', values[ka] - values[kb], units[ka], None, 'auto'))
                    n_auto += 1
    return cands, {'facts': len(raw), 'listed_derived': n_listed, 'auto_diff': n_auto, 'errors': errors}


class Index:
    """candidate values converted to one printed unit, sorted for range search"""
    def __init__(self, cands):
        self.cands = cands
        self.by_unit = {}

    def get(self, up):
        if up in self.by_unit:
            return self.by_unit[up]
        rows = []
        for i, c in enumerate(self.cands):
            if c.value is None:
                continue
            if up is None:
                rows.append((abs(c.value), i, ''))
                continue
            if c.unit is None:
                if c.kind == 'string':
                    rows.append((abs(c.value), i, 'unit?'))
                continue
            v, lab = convert(c.value, c.unit, up)
            if v is None:
                continue
            rows.append((abs(v), i, lab))
        rows.sort()
        self.by_unit[up] = ([r[0] for r in rows], rows)
        return self.by_unit[up]

    def search(self, x, tol, up):
        keys, rows = self.get(up)
        lo = bisect.bisect_left(keys, x - tol)
        hi = bisect.bisect_right(keys, x + tol)
        return rows[lo:hi]


# ---------------------------------------------------------------------------------------------------------------
# literature notes
# ---------------------------------------------------------------------------------------------------------------
def load_lit(lit_dir):
    pool_key, pool_all = defaultdict(list), []
    files = []
    if lit_dir and Path(lit_dir).is_dir():
        files = sorted([p for p in Path(lit_dir).glob('*.md') if p.name == 'REVIEW.md' or p.name.startswith('notes_')])
    key_re = re.compile(r'\b[a-z][a-z\-]*\d{4}[a-z0-9]+\b')
    num_re = re.compile(r'(?<![\w.])[−\-+]?\d+(?:\.\d+)?(?![\w]*\d)')
    for f in files:
        lines = f.read_text(encoding='utf-8', errors='replace').split('\n')
        head_keys = []
        block, b0 = [], 0

        def flush(block, b0, head_keys):
            text = '\n'.join(block)
            keys = set(key_re.findall(text)) | set(head_keys)
            for m in num_re.finditer(text):
                s = m.group().replace('−', '-').lstrip('+')
                d = len(s.split('.')[1]) if '.' in s else 0
                ln = b0 + text[:m.start()].count('\n') + 1
                rec = (abs(float(s)), d, f'{f.name}:L{ln}')
                pool_all.append(rec)
                for k in keys:
                    pool_key[k].append(rec)
        for i, ln in enumerate(lines):
            if re.match(r'^#{1,6}\s', ln):
                if block:
                    flush(block, b0, head_keys)
                block = []
                head_keys = key_re.findall(ln)
                continue
            if not ln.strip() or ln.lstrip().startswith('|'):
                if block:
                    flush(block, b0, head_keys)
                block = []
                if ln.lstrip().startswith('|'):
                    flush([ln], i, head_keys)
                continue
            if not block:
                b0 = i
            block.append(ln)
        if block:
            flush(block, b0, head_keys)
    return pool_key, pool_all, files


# ---------------------------------------------------------------------------------------------------------------
# manuscript numbers
# ---------------------------------------------------------------------------------------------------------------
class Num:
    __slots__ = ('text', 'value', 'decimals', 'neg', 'unit', 'a', 'b', 'line', 'approx', 'block', 'cites', 'ctx',
                 'status', 'match', 'alts', 'note', 'version', 'unit_src')

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))


APPROX_RE = re.compile(r'(?:about|approximately|approx\.|ca\.|nearly|roughly|around|some|almost|close to|over|under|'
                       r'below|above|less than|more than|up to|at least|at most|of order|≈|~|∼|<|>|≤|≥)\s*$', re.I)
SKIP_REF_RE = re.compile(
    r'\b(?:Figs?\.|Figures?|Tables?|Schemes?)(?:[   ]+|(?=\d))S?\d+[a-z]?(?:\s*(?:,|and|&|–|-|to)\s*S?\d+[a-z]?)*'
    r'|\b(?:Eqs?\.|Equations?|Eqns?\.)[   ]*\(?\d+\)?(?:\s*(?:,|and|&|–|-|to)\s*\(?\d+\)?)*'
    r'|\b(?:Sections?|Sects?\.|§)[   ]*\d+(?:\.\d+)*(?:\s*(?:,|and|&|–|-|to)\s*\d+(?:\.\d+)*)*'
    r'|\b(?:Refs?\.|references?)[   ]*\[?\d+\]?(?:\s*(?:,|and|–|-)\s*\d+)*')
METHOD_RE = re.compile(r'\b\d-\d+\+*G\**(?:\([^)]*\))?|\b(?:def2|cc|aug-cc|jun-cc|ma-def2)-\S+|\b10\.\d{4,}/\S+')
TOKEN_RE = re.compile(r'\d+(?:[.,]\d+)*')
SUP_DIGITS = str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺', '0123456789-+')
LOCANT_WORDS = re.compile(r'(?:amino|nitro|methyl|ethyl|chloro|bromo|fluoro|iodo|hydroxy|oxo|oxy|yl|ium|ide|ate|'
                          r'carbox|phen|benz|imidazol|triazol|pyrid|picr|nitrophen|trinitro|dinitro|di|tri|tetra|H\b|'
                          r'sulf|anilin)',   # SAPH: 4-sulfamoylanilinium, 4-sulfamoylphenyl
                          re.I)


def _sym(s):
    return re.sub(r'[\s*~^_{}\\]|\([^()]*\)', '', s)


def table_units(doc, blk, caption_text):
    """-> ({column index: unit}, default unit). Units come from the header cell ("D–H (Å)", "E / eV"), else from
    a caption exception naming that column's symbol ("(eV; *S* in eV⁻¹, *μ*~D~ in D)"), else from the caption's
    single default unit. A caption listing several units without saying which column is which gives no default."""
    lines = doc.text[blk.start:blk.end].split('\n')
    header = [c for c in lines[0].strip().strip('|').split('|')] if lines else []
    default, special = None, {}
    if caption_text:
        m = re.search(r'\(([^()]*(?:\([^()]*\)[^()]*)*)\)', caption_text)
        if m:
            plain = set()
            for part in re.split(r'[;,]', m.group(1)):
                mm = re.match(r'\s*(.+?)\s+in\s+(.+?)\s*$', part)
                if mm:
                    u = unit_in(mm.group(2))
                    if u:
                        special[_sym(mm.group(1))] = u
                    continue
                u = unit_in(part)
                if u:
                    plain.add(u)
            if len(plain) == 1:
                default = plain.pop()
    cols = {}
    for i, c in enumerate(header):
        u = unit_in(c) if re.search(r'\(|/|\[', c) else None
        if u is None and _sym(c) in special:
            u = special[_sym(c)]
        if u is None and default is not None and i > 0:
            u = default
        if u:
            cols[i] = u
    return cols, None


def unit_in(s):
    s = re.sub(r'[*~]', '', s)
    m = re.search(r'\(([^()]*)\)|/\s*(\S+)\s*$|\[([^\]]*)\]', s)
    cand = [g for g in (m.groups() if m else ()) if g] if m else []
    for c in cand + [s]:
        c = c.strip()
        for rx, u in PRINTED_UNITS:
            mm = rx.match(' ' + c)
            if mm and mm.end() == len(c) + 1:
                return u
    return None


def extract(doc):
    """-> (numbers, skipped Counter)"""
    from collections import Counter
    skipped = Counter()
    nums = []
    work, text = doc.work, doc.text
    # inline math: decimals only
    math_decimals = []
    i = 0
    while i < doc.n:
        if doc.cat[i] == MATH:
            j = i
            while j < doc.n and doc.cat[j] == MATH:
                j += 1
            seg = text[i:j]
            if seg.startswith('$$'):
                skipped['display math'] += len(re.findall(r'\d+', seg)) and 1
            else:
                for m in re.finditer(r'(?<![\w.^_{])([−\-]?)(\d+\.\d+)', seg):
                    math_decimals.append((i + m.start(2), i + m.end(2), m.group(1)))
                skipped['integers in math'] += len(re.findall(r'(?<![\d.])\d+(?![\d.])', seg))
            i = j
        else:
            i += 1
    caption_of = {}
    last_cap = None
    for b in doc.blocks:
        if b.kind == 'tcaption':
            last_cap = b
        if b.kind == 'table':
            cap = last_cap if last_cap and 0 < b.start - last_cap.end < 5 else None
            if cap is None:
                nxt = [x for x in doc.blocks if x.kind == 'tcaption' and 0 < x.start - b.end < 5]
                cap = nxt[0] if nxt else None
            caption_of[b.start] = cap
    for b in doc.blocks:
        if b.refs or b.kind in ('code', 'math', 'meta', 'refsdiv'):
            continue
        seg = work[b.start:b.end]
        skip_spans = [(m.start(), m.end(), 'figure/table/eq/section reference') for m in SKIP_REF_RE.finditer(seg)]
        skip_spans += [(m.start(), m.end(), 'method name / DOI') for m in METHOD_RE.finditer(seg)]
        if b.kind == 'heading':
            m = re.match(r'#+\s+\d+(?:\.\d+)*\.?', seg)
            if m:
                skip_spans.append((0, m.end(), 'section number'))
        cols, default = ({}, None)
        cap = None
        if b.kind == 'table':
            cap = caption_of.get(b.start)
            cols, default = table_units(doc, b, doc.text[cap.start:cap.end] if cap else '')
        cites_block = [c for c in doc.citations if c.lit_keys and b.start <= c.start < b.end]
        sents = doc.sentences(b) if b.kind in ('abstract', 'para', 'list', 'quote', 'figure', 'tcaption', 'title') else []
        extra = [(a - b.start, e - b.start, s) for a, e, s in math_decimals if b.start <= a < b.end]
        found = []
        for m in TOKEN_RE.finditer(seg):
            found.append((m.start(), m.end(), None))
        for a, e, s in extra:
            found.append((a, e, s))
        found.sort()
        block_nums = []
        consumed = -1                              # end of the last number (scientific notation, unit)
        for s0, s1, math_sign in found:
            if s0 < consumed:
                continue
            tok = seg[s0:s1] if math_sign is None else doc.text[b.start + s0:b.start + s1]
            A, B = b.start + s0, b.start + s1
            prev = seg[s0 - 1] if s0 > 0 else ' '
            prev2 = seg[s0 - 2] if s0 > 1 else ' '
            nxt = seg[s1] if s1 < len(seg) else ' '
            nxt2 = seg[s1 + 1] if s1 + 1 < len(seg) else ' '
            why = None
            for x0, x1, reason in skip_spans:
                if x0 <= s0 < x1:
                    why = reason
                    break
            if why is None and math_sign is None:
                if prev.isalpha() or prev in '_^~\'′' or (prev == '.' and prev2.isdigit()):
                    why = 'atom label / formula / subscript'
                elif prev in '-‐‑+−' and prev2.isalnum():
                    why = 'compound token (P-1, 6-311, -3-ium, HOMO−1, LUMO+1)'
                elif prev in ')]' or (prev in '([' and prev2.isdigit()):
                    why = 'standard uncertainty / charge'
                elif nxt == ':' and nxt2.isdigit() or (prev == ':' and prev2.isdigit()):
                    why = 'ratio'
                elif ',' in tok and not re.fullmatch(r'\d{1,3}(?:,\d{3})+(?:\.\d+)?', tok):
                    why = 'locant list'
                elif nxt == 'H' and nxt2 in '-,‐':
                    why = 'locant'
                elif nxt in '-‐‑' and nxt2.isalpha() and '.' not in tok and int(tok) < 10 \
                        and LOCANT_WORDS.match(seg[s1 + 1:s1 + 12]):
                    why = 'locant'
                elif re.match(r'(?:st|nd|rd|th)\b', seg[s1:s1 + 3]):
                    why = 'ordinal'
                elif nxt.isalpha() and (('.' not in tok and re.match(r'D(?![\w])', seg[s1:]))
                                        or not re.match(r'(?:Å|eV|kcal|kJ|cm|nm|pm|mm|m?Ha|K|D|a\.u\.|au\b|Da|ppm)',
                                                        seg[s1:])):
                    why = 'label (3D, 2a, 4c)'
            if why is None and math_sign is None and (
                    re.search(r'(?<![\w])\*?[xyz]\*?\s*[+−-]\s*$', seg[max(0, s0 - 6):s0]) or
                    re.match(r'\s*[+−-]\s*\*?[xyz]\*?(?![\w])', seg[s1:s1 + 6])):
                why = 'symmetry operator'
            if why is None and re.fullmatch(r'(19[5-9]\d|20[0-3]\d)', tok):
                keyyear = any(tok in k for c in doc.citations if c.lit_keys and abs(c.start - A) < 400 for k in c.lit_keys)
                if re.search(r'et al\.?,?\s*\(?$', seg[max(0, s0 - 12):s0]) or keyyear:
                    why = 'year in citation'
            if why:
                skipped[why] += 1
                continue
            version = tok.count('.') > 1
            neg = False
            if math_sign:
                neg = math_sign in '−-'
            elif prev in '−-' and (s0 < 2 or prev2 in ' \t\n([=≈~<>≤≥±/;:,  '):
                neg = True
            clean = tok.replace(',', '')
            dec = len(clean.split('.')[1]) if ('.' in clean and not version) else 0
            try:
                val = float(clean) if not version else None
            except ValueError:
                val = None
            # scientific notation
            rest = seg[s1:s1 + 20]
            sm = re.match(r'\s*[×x]\s*10(?:\^\{?([−\-+]?\d+)\}?\^?|([⁻⁺]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+))', rest) or \
                re.match(r'[eE]([−\-+]?\d+)\b', rest)
            end = s1
            pw = re.match(r'\^\{?([−\-+]?\d+)\}?\^?|([⁻⁺]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+)', rest) if tok == '10' else None
            if pw and not sm:                              # a bare power of ten: 10⁻³, 10^-3^
                exp = int((pw.group(1) or pw.group(2).translate(SUP_DIGITS)).replace('−', '-'))
                val, dec, end = 10.0 ** exp, -exp, s1 + pw.end()
            if sm and val is not None:
                exp = sm.group(1) if sm.group(1) else sm.group(2).translate(SUP_DIGITS)
                exp = int(exp.replace('−', '-'))
                val = val * 10 ** exp
                dec = dec - exp
                end = s1 + sm.end()
            end_num = end
            unit = None
            um = None
            if not version:
                after = seg[end:end + 25]
                for rx, u in PRINTED_UNITS:
                    um = rx.match(after)
                    if um:
                        unit = u
                        end += um.end()
                        break
            unit_src = 'text' if unit else None
            if unit is None and b.kind == 'table':
                line_start = seg.rfind('\n', 0, s0) + 1
                col = seg[line_start:s0].count('|') - 1
                unit = cols.get(col, default)
                unit_src = 'table' if unit else None
            approx = bool(APPROX_RE.search(seg[max(0, s0 - 20):s0 - (1 if neg else 0)]))
            # citations governing this number
            if b.kind == 'table':
                ls = seg.rfind('\n', 0, s0) + 1
                le = seg.find('\n', s0)
                le = len(seg) if le < 0 else le
                cites = [k for c in cites_block if b.start + ls <= c.start < b.start + le for k in c.lit_keys]
                if cap is not None:
                    cites += [k for c in doc.citations if c.lit_keys and cap.start <= c.start < cap.end for k in c.lit_keys]
            else:
                sa, se = next(((x, y) for x, y in sents if x <= A < y + 1), (b.start, b.end))
                cites = [k for c in cites_block if sa <= c.start <= se + 2 for k in c.lit_keys]
            consumed = end
            shown = re.sub(r'\s+', ' ', doc.text[A:b.start + end_num]).strip()
            block_nums.append(Num(text=('−' if neg else '') + shown, value=val, decimals=dec, neg=neg, unit=unit,
                                  a=A - (1 if neg else 0), b=b.start + end, line=doc.line_of(A), approx=approx,
                                  block=b.kind, cites=cites, ctx=doc.context(A - (1 if neg else 0), b.start + end, 40),
                                  version=version, unit_src=unit_src))
        # unit inheritance across ranges and lists: "1.052 and 1.702 Å", "−0.07 to +0.09 a.u.", "0.03 ± 0.01 Å"
        for k in range(len(block_nums) - 2, -1, -1):
            x, y = block_nums[k], block_nums[k + 1]
            if x.unit is None and y.unit is not None:
                gap = work[x.b:y.a]
                if re.fullmatch(r'[\s  ]*(?:–|-|to|and|or|,|, and|/|±|and\s+[+−-]?|,\s*[+−-]?|to\s+[+−-]?)?'
                                r'[\s  ]*[+]?', gap):
                    x.unit, x.unit_src = y.unit, 'inherited'
        nums.extend(block_nums)
    return nums, skipped


def tolerance(n):
    tol = 0.5 * 10 ** (-n.decimals)
    if n.approx and n.decimals == 0 and n.value:
        t = len(n.text) - len(n.text.rstrip('0'))
        if t:
            tol = max(tol, 0.5 * 10 ** t)
    return tol * (1 + 1e-9) + 1e-12


CTX_CASE = {
    'gas': re.compile(r'gas[- ]phase|\bgas\b|vacuo|vacuum|isolated', re.I),
    'pcm': re.compile(r'pcm|water|continuum|solvat|solvent|aqueous', re.I),
    'xray': re.compile(r'x-ray|crystal|\bcif\b|diffraction|experiment', re.I),
    'periodic': re.compile(r'periodic|plane[- ]wave|quantum espresso|\bPBE-D3\b|unit cell', re.I),
    'saddle': re.compile(r'saddle|transition', re.I),
}
GREEK = {'η': 'hardness eta', 'ω': 'electrophilicity omega', 'χ': 'electronegativity chi', 'μ': 'mu dipole potential',
         'ρ': 'rho density', 's': 'softness', 'δ': 'delta', 'λ': 'lambda'}


def score(n, c, lab, ctx_words, ctx_cases, y=None, tol=None):
    s = {'value': 3.0, 'member': 3.0, 'string': 2.5, 'derived': 2.0, 'auto': 1.0}[c.kind]
    if y is not None and tol:
        s += 1.0 - min(1.0, abs(abs(n.value) - y) / tol)       # the closer to the printed value, the better
    if n.unit is not None and c.unit == n.unit:
        s += 2
    elif lab and lab != 'unit?':
        s += 1.5
    elif lab == 'unit?':
        s -= 0.5
    if c.case and c.case in ctx_cases:
        s += 1
    toks = [t for t in re.split(r'[._\[\]\s=−-]+', c.key.lower()) if len(t) >= 3]
    s += min(2.0, 0.5 * sum(1 for t in toks if t in ctx_words))
    if c.value is not None and (c.value < 0) == n.neg:
        s += 0.3
    return s


def audit(man, facts_path, lit_dir, auto_diff=True):
    doc = Doc.from_file(man)
    cands, finfo = load_facts(facts_path, auto_diff)
    idx = Index(cands)
    pool_key, pool_all, lit_files = load_lit(lit_dir)
    nums, skipped = extract(doc)
    versions = {c.text for c in cands if c.kind == 'version'}
    for n in nums:
        ctx_text = n.ctx.lower()
        for g, words in GREEK.items():
            if g in n.ctx:
                ctx_text += ' ' + words
        ctx_words = set(re.findall(r'[a-z]{3,}', ctx_text))
        ctx_cases = {k for k, rx in CTX_CASE.items() if rx.search(n.ctx)}
        if n.version:
            n.status = 'matched' if n.text in versions else 'unmatched'
            n.match = ('version string', '')
            if n.status == 'unmatched':
                n.status = None
        elif n.value is not None:
            x = abs(n.value)
            tol = tolerance(n)
            rows = idx.search(x, tol, n.unit)
            if rows:
                best = max(rows, key=lambda r: score(n, cands[r[1]], r[2], ctx_words, ctx_cases, r[0], tol))
                c = cands[best[1]]
                n.status = 'matched'
                n.alts = len({cands[r[1]].key for r in rows}) - 1
                notes = []
                if best[2] and best[2] != 'unit?':
                    notes.append(best[2])
                if best[2] == 'unit?':
                    notes.append('fact has no unit')
                if c.value is not None and (c.value < 0) != n.neg:
                    notes.append('|value|' if c.value < 0 else 'sign flipped')
                if c.kind in ('derived', 'auto'):
                    notes.append('derived' if c.kind == 'derived' else 'auto difference')
                if n.approx and n.decimals == 0 and tol > 0.5000001:
                    notes.append(f'approx. ±{tol:.3g}')
                sig = len(re.sub(r'^[0.]+', '', re.sub(r'[^0-9.]', '', n.text.split('×')[0])).replace('.', ''))
                if sig <= 1 and n.alts >= 2:
                    notes.append(f'weak: {sig} significant figure, {n.alts + 1} candidates')
                vdisp = f'{c.value:.6g}' if c.kind != 'string' else c.text
                n.match = (c.key, f'{vdisp} {c.unit or ""}'.strip())
                n.note = ', '.join(notes)
        if n.status is None and n.cites:
            x = abs(n.value) if n.value is not None else None
            tol = tolerance(n) if x is not None else 0
            near = [r for k in n.cites for r in pool_key.get(k, []) if x is not None and abs(r[0] - x) <= tol]
            if near:
                n.status, n.match, n.note = 'literature', (', '.join(sorted(set(n.cites))), near[0][2]), 'beside key'
            else:
                anyw = [r for r in pool_all if x is not None and abs(r[0] - x) <= tol]
                if anyw:
                    n.status, n.match, n.note = 'literature', (', '.join(sorted(set(n.cites))), anyw[0][2]), \
                        'elsewhere in notes'
        if n.status is None:
            n.status = 'unmatched'
    return doc, nums, skipped, finfo, lit_files


def resolve(arg, a):
    if arg in PAPER_IDS:
        base = PAPERS / arg
    else:
        base = Path(arg)
        if base.is_file():
            base = base.parent
    man = Path(a.manuscript) if a.manuscript else base / 'manuscript.md'
    # The facts dossier is built by paper/main/data/build_facts.py into paper/main/data/facts.json (one file for
    # the paper and the SI); paper/data/facts.json is the 2-entry seed of the skeleton and is no longer used.
    facts = Path(a.facts) if a.facts else PAPERS / 'main' / 'data' / 'facts.json'
    lit = Path(a.lit) if a.lit else LIT
    return man, facts, lit


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('paper', help='main, si, or a directory holding manuscript.md')
    ap.add_argument('--manuscript')
    ap.add_argument('--facts')
    ap.add_argument('--lit')
    ap.add_argument('--no-auto-diff', action='store_true', help='do not add gas/pcm/xray counterpart differences')
    ap.add_argument('--quiet', action='store_true', help='print UNMATCHED and the summary only')
    ap.add_argument('--json', help='write the result as JSON')
    a = ap.parse_args(argv)
    man, facts, lit = resolve(a.paper, a)
    for p in (man, facts):
        if not p.exists():
            print(f'numbers_audit: {p} not found', file=sys.stderr)
            return 2
    doc, nums, skipped, finfo, lit_files = audit(man, facts, lit, not a.no_auto_diff)
    print(f'== numbers_audit: {rel(man)}\n   facts: {rel(facts)} ({finfo["facts"]} entries, {finfo["listed_derived"]} '
          f'listed derived, {finfo["auto_diff"]} automatic differences); literature notes: '
          f'{", ".join(p.name for p in lit_files) or "none"}')
    for e in finfo['errors']:
        print(f'   WARNING {e}')
    groups = defaultdict(list)
    for n in nums:
        groups[n.status].append(n)

    def fmt_num(n):
        u = f' {n.unit}' if n.unit else ''
        src = {'inherited': '*', 'table': '†'}.get(n.unit_src, '')
        return f'{n.text}{u}{src}'
    if not a.quiet:
        print(f'\n-- MATCHED ({len(groups["matched"])})   [* unit inherited from the next number, † unit from table header/caption]')
        for n in groups['matched']:
            alt = f' (+{n.alts} other)' if n.alts else ''
            note = f'  [{n.note}]' if n.note else ''
            print(f'  L{n.line:<5d} {fmt_num(n):22s} → {n.match[0]} = {n.match[1]}{note}{alt}')
        print(f'\n-- LITERATURE ({len(groups["literature"])})')
        for n in groups['literature']:
            print(f'  L{n.line:<5d} {fmt_num(n):22s} [@{n.match[0]}] → {n.match[1]} ({n.note})')
    print(f'\n-- UNMATCHED ({len(groups["unmatched"])})')
    for n in groups['unmatched']:
        cite = f' [cites @{", @".join(n.cites)}]' if n.cites else ''
        print(f'  L{n.line:<5d} {fmt_num(n):22s} {n.ctx}{cite}')
    sk = ', '.join(f'{k} {v}' for k, v in sorted(skipped.items(), key=lambda kv: -kv[1]) if v)
    print(f'\n-- skipped (not audited): {sum(skipped.values())} ({sk or "none"})')
    nun = len(groups['unmatched'])
    print(f'-- summary: {len(nums)} numbers: {len(groups["matched"])} matched, {len(groups["literature"])} literature, '
          f'{nun} unmatched -> {"FAIL" if nun else "PASS"}')
    if a.json:
        out = [{k: getattr(n, k) for k in ('line', 'text', 'value', 'unit', 'status', 'match', 'note', 'cites', 'ctx')}
               for n in nums]
        Path(a.json).write_text(json.dumps({'numbers': out, 'skipped': dict(skipped)}, indent=1, ensure_ascii=False),
                                encoding='utf-8')
    return 1 if nun else 0


if __name__ == '__main__':
    sys.exit(main())
