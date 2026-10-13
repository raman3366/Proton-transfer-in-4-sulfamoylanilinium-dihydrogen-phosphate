#!/usr/bin/env python3
"""lint_prose.py -- check a manuscript against paper/_shared/STYLE.md (standard library only).

Usage
  python3 paper/_shared/tools/lint_prose.py <manuscript.md | main | si | paper dir> [...]
          [--style paper/_shared/STYLE.md] [--json report.json] [--quiet] [--bigrams 30]

  main / si are shorthand for paper/main/manuscript.md and paper/si/manuscript.md (adapted from cqc-dft).

The banned words and phrases are parsed from the STYLE.md section "Words and constructions to avoid" (the
list lives only there). Items are matched case-insensitively on word boundaries, with simple inflections
(-s/-es/-d/-ed/-ing, e-drop, y->ies, consonant doubling; the first word of verb phrases such as "plays a ...
role", "sheds light", "serves as"; an optional plural on the last word of a phrase). Qualifiers written in
STYLE.md are honoured: "(max once per paper)", "(at most once)", "(max twice)" set a limit per manuscript;
"(as a verb)" only flags verb use; "(unless statistical)" / "(unless a statistical test was done)" become a
warning when the sentence mentions a statistical test. "electron-rich" and similar compounds are not "rich".

Checks (HARD = exit status 1)
  HARD  banned        banned word/phrase, or an allowed-N-times phrase used more than N times
  HARD  em-dash       — or ― (or pandoc '---') outside code and math
  HARD  schiff        the word Schiff anywhere in the text
  HARD  symmetry      internal symmetry suffixes (#1_545, N4#2, 1_565) outside comments/code
  HARD  floats        figure/table defined but never referenced; reference to an undefined figure/table
  warn  floats        figures/tables first mentioned out of order; unlabelled captioned figure
  warn  openers       consecutive sentences with the same first two words; a two-word opener used > 3 times
                      in one section; consecutive paragraphs with the same opener
  warn  restatement   paragraph-final "Thus, the / This indicates that / These results / In summary / Overall ..."
  warn  citations     citation bracket with > 3 keys (also two adjacent brackets together > 3)
  warn  rhythm        sentence-length coefficient of variation < 0.35 in a section (>= 5 sentences)
  warn  no-number     paragraph in Results and discussion without a single digit
  warn  units         kcal/mol, cm-1, cm^-1, mol-1, Å2, 'A' or 'angstrom' for Å, number glued to its unit,
                      'deg' for °, 'Debye', mixed 'a.u.'/'au'
  warn  style         !, ?, **bold**, bullet lists in the body, stacked hedges, "on the one hand", '...' or '⋯'
                      instead of ···, hyphen instead of en dash in bonds (N4-H4), abstract 180-250 words,
                      conclusions 250-400 words, 5-6 keywords
  warn  docking       SAPH docking wording: "binding energy/affinity", K_i/K_d/IC50, "strong binding",
                      "potent inhibitor", "drug candidate", "lead", and SAPH / the salt / the ion pair "docked"
  info  bigrams       the most frequent content bigrams (repetition check)
"""
import argparse
import json
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mdtext import (Doc, SHARED, XREF, CITE, ATTR, COMMENT, HTML, CODE, META, WORD_RE,  # noqa: E402
                    resolve_manuscript, rel, strip_number)

STYLE = SHARED / 'STYLE.md'
HARD_CHECKS = ('banned', 'em-dash', 'schiff', 'symmetry', 'floats-hard')

STATS_RE = re.compile(r'statistic|p[- ]value|\bp\s*[<≤=]|t-test|confidence interval|standard deviation|anova|'
                      r'regression|correlation coefficient|\br(?:\^?2\^?|²)\b|significance level|chi-squared', re.I)
VERB_BEFORE = re.compile(r'\b(to|we|can|could|will|would|may|might|should|must|they|which|that|and|or|also)\s+$', re.I)
VERB_AFTER = re.compile(r'^\s+(the|a|an|this|these|those|our|its|their|it|them|such)\b', re.I)
NO_INFLECT_FIRST = {'this', 'has', 'was', 'its', 'his', 'does', 'is', 'as', 'us', 'thus', 'always', 'less'}

STOP = set('''a about above after again against all almost along also although always am among an and another any
are around as at be because been before being below between both but by can cannot could did do does doing down
during each either else even ever every for from further had has have having he her here hers him his how however i
if in into is it its itself just least less let like made make many may me might more most much must my near neither
no nor not now of off on once one only or other others otherwise our ours out over own per rather same she should
since so some such than that the their theirs them then there therefore these they this those though through thus
to too toward towards under until up upon us very via was we well were what whatever when where whereas whether which
while who whom whose why will with within without would yet you your fig figs table tables eq eqs section sections
see i.e e.g ca vs both whereas here'''.split())


# --------------------------------------------------------------------------------------------------------------
# STYLE.md banned list
# --------------------------------------------------------------------------------------------------------------
class Rule:
    def __init__(self, label, regex, limit=None, mode='ban', note=''):
        self.label, self.regex, self.limit, self.mode, self.note = label, regex, limit, mode, note

    def __repr__(self):
        return f'Rule({self.label!r}, limit={self.limit}, mode={self.mode})'


def _split_top(s):
    items, buf, depth, inq = [], [], 0, False
    for ch in s:
        if ch in '"“”':
            inq = not inq if ch == '"' else (ch == '“')
            buf.append(ch)
            continue
        if not inq:
            if ch == '(':
                depth += 1
            elif ch == ')':
                depth = max(0, depth - 1)
            elif ch == ',' and depth == 0:
                items.append(''.join(buf))
                buf = []
                continue
        buf.append(ch)
    if ''.join(buf).strip():
        items.append(''.join(buf))
    return [i.strip() for i in items if i.strip()]


def _inflect(w):
    out = {w, w + 's', w + 'es', w + 'd', w + 'ed', w + 'ing'}
    if w.endswith('e'):
        out |= {w[:-1] + 'ing', w[:-1] + 'ed'}
    if w.endswith('y'):
        out |= {w[:-1] + 'ies', w[:-1] + 'ied'}
    if re.search(r'[^aeiou][aeiou][bdglmnprt]$', w):
        out |= {w + w[-1] + 'ed', w + w[-1] + 'ing'}
    return out


def _word_alts(word, inflect):
    """one word of an item -> list of literal alternatives"""
    alts = set()
    for part in word.split('/'):
        m = re.fullmatch(r'([\w\'-]+)\((\w+)\)', part)
        if m:
            base, suf = m.group(1), m.group(2)
            forms = {base, base + suf}
            if inflect:
                forms |= _inflect(base)
        else:
            forms = {part}
            if inflect:
                forms |= _inflect(part)
        alts |= forms
    return sorted(alts, key=len, reverse=True)


def _alt_re(alts):
    out = []
    for a in alts:
        e = re.escape(a).replace(r'\-', '-')
        e = e.replace('-', r'[-‐‑\s]+')
        out.append(e)
    return '(?:' + '|'.join(out) + ')'


def parse_style(path=STYLE):
    text = Path(path).read_text(encoding='utf-8')
    m = re.search(r'^##\s+Words and constructions to avoid\s*$(.*?)(?=^##\s|\Z)', text, re.M | re.S)
    if not m:
        raise SystemExit(f'{path}: section "Words and constructions to avoid" not found')
    body = re.split(r'\n\s*Also avoid', m.group(1))[0]
    colon = re.search(r':\s*\n', body)
    lst = body[colon.end():] if colon else body
    lst = ' '.join(lst.split())
    rules, seen = [], set()
    for raw in _split_top(lst):
        s = raw.strip().rstrip('.')
        quals = []
        while True:
            q = re.search(r'\s+\(([^()]*)\)\s*$', s)
            if not q:
                break
            quals.insert(0, q.group(1).strip())
            s = s[:q.start()].rstrip()
        s = s.strip().strip('"“”').strip()
        if not s:
            continue
        key = s.lower()
        if key in seen:
            continue
        seen.add(key)
        limit, mode, note = None, 'ban', '; '.join(quals)
        for q in quals:
            ql = q.lower()
            mm = re.search(r'(?:max|at most)\s+(once|twice|thrice|\d+)', ql)
            if mm:
                limit = {'once': 1, 'twice': 2, 'thrice': 3}.get(mm.group(1)) or int(mm.group(1))
            elif ql.startswith('unless'):
                mode = 'unless_stat'
            elif 'as a verb' in ql:
                mode = 'verb'
        words = [w.rstrip(',') for w in s.split()]
        parts = []
        for i, w in enumerate(words):
            if w in ('...', '…'):
                parts.append(None)
                continue
            if len(words) == 1:
                alts = _word_alts(w, True)
            elif i == 0 and len(w) > 3 and w.lower().endswith('s') and w.lower() not in NO_INFLECT_FIRST \
                    and not w[0].isupper():
                base = w[:-1]
                alts = sorted({w} | _inflect(base), key=len, reverse=True)
            elif i == len(words) - 1 and re.fullmatch(r'[a-z]{3,}', w):
                alts = [w + 'es', w + 's', w]
            else:
                alts = _word_alts(w, False)
            parts.append(_alt_re(alts))
        rx = ''
        for i, p in enumerate(parts):
            if p is None:
                rx += r'[^.;:!?]{1,120}?'
                continue
            if i > 0 and parts[i - 1] is not None:
                rx += r',?\s+'
            elif i > 0:
                rx += r'\b'
            rx += p
        regex = re.compile(r'(?<![\w\-])' + rx + r'(?![\w])', re.I)
        rules.append(Rule(s, regex, limit, mode, note))
    return rules


# --------------------------------------------------------------------------------------------------------------
# findings
# --------------------------------------------------------------------------------------------------------------
class Findings:
    def __init__(self, doc):
        self.doc = doc
        self.items = []          # (severity, check, line, message, context)

    def add(self, sev, check, a, b, msg, ctx=True):
        msg = re.sub(r'\s*\n\s*', ' ', msg)
        line = self.doc.line_of(a) if a is not None else 0
        c = self.doc.context(a, b) if (ctx and a is not None) else ''
        self.items.append((sev, check, line, msg, c))

    def hard(self):
        return [i for i in self.items if i[0] == 'HARD']


def body_blocks(doc, kinds=('title', 'abstract', 'keywords', 'heading', 'para', 'figure', 'tcaption', 'table',
                            'list', 'quote')):
    return [b for b in doc.blocks if b.kind in kinds and not b.refs]


def in_blocks(blocks, off):
    for b in blocks:
        if b.start <= off < b.end:
            return b
    return None


_SENT_CACHE = {}


def sentence_at(doc, off):
    """-> (start, end) of the sentence containing offset off (the block if it has no sentences)"""
    b = doc.block_at(off)
    if b is None:
        return off, off
    key = (id(doc), b.start)
    if key not in _SENT_CACHE:
        _SENT_CACHE[key] = doc.sentences(b) or [(b.start, b.end)]
    for x, y in _SENT_CACHE[key]:
        if x <= off <= y:
            return x, y
    return b.start, b.end


def check_banned(doc, F, rules):
    blocks = body_blocks(doc)
    counts = defaultdict(list)
    raw_hits = []
    for r in rules:
        for m in r.regex.finditer(doc.work):
            if in_blocks(blocks, m.start()) is not None:
                raw_hits.append((m.start(), m.end(), r, m))
    # a hit lying inside a longer hit of another rule ("pivotal" inside "plays a pivotal role") is dropped
    hits = [h for h in raw_hits if not any(o is not h and o[2] is not h[2] and o[0] <= h[0] and h[1] <= o[1]
                                           and (o[1] - o[0]) > (h[1] - h[0]) for o in raw_hits)]
    for a, e, r, m in hits:
        low = m.group().lower()
        if low.startswith('rich') and a > 0 and doc.work[a - 1] in '-‐‑':
            continue
        if r.mode == 'verb':
            before = doc.work[max(0, a - 20):a]
            after = doc.work[e:e + 12]
            if low == r.label.lower() and not (VERB_BEFORE.search(before) or VERB_AFTER.match(after)):
                F.add('warn', 'banned', a, e, f'"{m.group()}" ({r.note}): check it is not a verb')
                continue
        if r.mode == 'unless_stat':
            x, y = sentence_at(doc, a)
            if STATS_RE.search(doc.work[x:y]):
                F.add('warn', 'banned', a, e, f'"{m.group()}" ({r.note}): statistical context found, check')
                continue
        counts[r.label].append((a, e, m.group(), r))
    for label, hits in counts.items():
        r = hits[0][3]
        hits.sort()
        if r.limit is not None:
            if len(hits) > r.limit:
                for a, e, g, _ in hits[r.limit:]:
                    F.add('HARD', 'banned', a, e, f'"{g}" used {len(hits)} times (STYLE allows {r.limit})')
            continue
        for a, e, g, _ in hits:
            note = f' ({r.note})' if r.note else ''
            F.add('HARD', 'banned', a, e, f'banned "{r.label}"{note}: "{g}"')


def check_dashes_schiff_symmetry(doc, F):
    blocks = body_blocks(doc)
    for m in re.finditer(r'[—―]', doc.work):
        if in_blocks(blocks, m.start()):
            F.add('HARD', 'em-dash', m.start(), m.end(), 'em dash: use commas, parentheses, a semicolon or a new sentence')
    for b in blocks:
        if b.kind == 'table':
            continue
        seg = doc.work[b.start:b.end]
        for m in re.finditer(r'(?<!-)---(?!-)', seg):
            F.add('HARD', 'em-dash', b.start + m.start(), b.start + m.end(), "'---' becomes an em dash in pandoc")
    for m in re.finditer(r'schiff', doc.work, re.I):
        if in_blocks(blocks, m.start()):
            F.add('HARD', 'schiff', m.start(), m.end(), 'never call these compounds Schiff bases')
    for m in re.finditer(r'schiff', doc.text, re.I):   # also in captions/alt text/attributes, not in comments
        if doc.cat[m.start()] in (ATTR, XREF, HTML) :
            F.add('HARD', 'schiff', m.start(), m.end(), 'Schiff in a label, path or attribute')
    pat = re.compile(r'(?<=[\w′\'])#\d+(?:_\d{3})?|#\d+_\d{3}\b|\b\d_\d{3}\b|\b[A-Za-z]{1,2}\d+[A-Za-z]?_\d_\d{3}\b')
    for m in pat.finditer(doc.text):
        if doc.cat[m.start()] in (COMMENT, CODE, META, ATTR, XREF):
            continue                           # comments, code, YAML keys, file paths and labels are never printed
        blk = doc.block_at(m.start())
        if blk is not None and blk.refs:
            continue
        F.add('HARD', 'symmetry', m.start(), m.end(),
              'internal symmetry suffix: describe the symmetry-generated fragment in words')


def check_floats(doc, F):
    defs = doc.float_defs()
    refs = doc.float_refs()
    names = {'fig': 'Fig.', 'tbl': 'Table'}
    for kind in ('fig', 'tbl'):
        items = defs[kind]
        lab2num = {lab: i + 1 for i, (lab, _) in enumerate(items) if lab}
        blocks = [b for _, b in items]
        first = {}
        for k, num, lab, off, literal in refs:
            if k != kind:
                continue
            if lab is not None:
                if lab not in lab2num:
                    F.add('HARD', 'floats-hard', off, off + len(lab) + 1, f'reference to undefined label @{lab}')
                    continue
                num = lab2num[lab]
            if num is None:
                continue
            if num < 1 or num > len(items):
                F.add('HARD', 'floats-hard', off, off + 6,
                      f'{names[kind]} {num} referenced but only {len(items)} defined')
                continue
            own = blocks[num - 1]
            if own.start <= off < own.end:
                continue                     # a caption mentioning its own number
            first.setdefault(num, off)
        for i, (lab, b) in enumerate(items):
            num = i + 1
            if lab is None and kind == 'fig':
                F.add('warn', 'floats', b.start, b.start + 20, f'Fig. {num} has no {{#fig:label}}')
            if num not in first:
                F.add('HARD', 'floats-hard', b.start, b.start + 30,
                      f'{names[kind]} {num} ({lab or "unlabelled"}) is never referenced in the text')
        for num in range(2, len(items) + 1):
            if num in first and (num - 1) in first and first[num] < first[num - 1]:
                F.add('warn', 'floats', first[num], first[num] + 8,
                      f'{names[kind]} {num} is first mentioned (L{doc.line_of(first[num])}) before '
                      f'{names[kind]} {num - 1} (L{doc.line_of(first[num - 1])})')
            elif num in first and (num - 1) not in first:
                pass


def sentence_units(doc):
    """-> {section: [(start, end, block)]} for running prose (abstract + paragraphs + lists), no back matter"""
    out = defaultdict(list)
    order = []
    for b in doc.blocks:
        if b.kind not in ('abstract', 'para', 'list', 'quote') or b.back or b.refs:
            continue
        sec = b.sec_fine or b.sec_top or '(no section)'
        if sec not in out:
            order.append(sec)
        for a, e in doc.sentences(b):
            out[sec].append((a, e, b))
    return order, out


def opener(doc, a, e):
    w = [x.lower() for x in WORD_RE.findall(doc.work[a:e])[:2]]
    return ' '.join(w) if len(w) == 2 else None


RESTATE_RE = re.compile(r'^[*_\s]*(thus,?\s+the|this\s+indicates\s+that|these\s+results|in\s+summary|overall|'
                        r'taken\s+together|collectively|in\s+short|together,?\s+these|these\s+findings|'
                        r'this\s+shows\s+that|this\s+confirms\s+that|hence,?\s+the|in\s+conclusion)\b', re.I)


def check_sentences(doc, F, stats_out):
    order, secs = sentence_units(doc)
    for sec in order:
        sents = secs[sec]
        ops = [opener(doc, a, e) for a, e, _ in sents]
        for i in range(1, len(sents)):
            if ops[i] and ops[i] == ops[i - 1]:
                a, e, _ = sents[i]
                F.add('warn', 'openers', a, min(e, a + 30), f'consecutive sentences open with "{ops[i]}"')
        cnt = Counter(o for o in ops if o)
        for o, n in cnt.items():
            if n > 3:
                lines = [doc.line_of(sents[i][0]) for i in range(len(sents)) if ops[i] == o]
                a = sents[ops.index(o)][0]
                F.add('warn', 'openers', a, a + len(o), f'opener "{o}" used {n} times in "{sec}" '
                      f'(lines {", ".join(map(str, lines))})', ctx=False)
        lens = [len(WORD_RE.findall(doc.work[a:e])) for a, e, _ in sents]
        if lens:
            mean = statistics.fmean(lens)
            sd = statistics.stdev(lens) if len(lens) > 1 else 0.0
            cv = sd / mean if mean else 0.0
            stats_out.append((sec, len(lens), mean, sd, cv))
            if len(lens) >= 5 and cv < 0.35:
                a = sents[0][0]
                F.add('warn', 'rhythm', a, a, f'sentence-length CV {cv:.2f} < 0.35 in "{sec}" '
                      f'(n={len(lens)}, mean {mean:.1f} words, sd {sd:.1f}): vary sentence length', ctx=False)
    # paragraph-final restatements and paragraph openers
    prev_op, prev_top = None, None
    for b in doc.blocks:
        if b.kind != 'para' or b.back or b.refs:
            if b.kind == 'heading':
                prev_op = None
            continue
        ss = doc.sentences(b)
        if not ss:
            continue
        a, e = ss[-1]
        if len(ss) > 1 and RESTATE_RE.match(doc.work[a:e]):
            F.add('warn', 'restatement', a, min(e, a + 40), 'paragraph ends with a summary sentence that restates it')
        op = opener(doc, *ss[0])
        if op and op == prev_op and b.sec_top == prev_top:
            F.add('warn', 'openers', ss[0][0], ss[0][0] + len(op), f'consecutive paragraphs open with "{op}"')
        prev_op, prev_top = op, b.sec_top


def check_citations(doc, F):
    br = [c for c in doc.citations if c.bracket and c.lit_keys]
    for c in br:
        if len(c.lit_keys) > 3:
            F.add('warn', 'citations', c.start, c.end, f'{len(c.lit_keys)} references in one bracket (max 3)')
    for c1, c2 in zip(br, br[1:]):
        gap = doc.text[c1.end:c2.start]
        if len(gap) <= 3 and re.fullmatch(r'[\s,;]*', gap) and len(c1.lit_keys) + len(c2.lit_keys) > 3:
            F.add('warn', 'citations', c1.start, c2.end,
                  f'adjacent brackets cite {len(c1.lit_keys) + len(c2.lit_keys)} references together')


def check_numbers_present(doc, F):
    for b in doc.blocks:
        if b.kind not in ('para', 'list') or b.back or b.refs:
            continue
        if not strip_number(b.sec_top).lower().startswith('results'):
            continue
        raw = ''.join(ch if doc.cat[i] not in (CITE, XREF, ATTR, COMMENT, HTML) else ' '
                      for i, ch in enumerate(doc.text[b.start:b.end], start=b.start))
        if not re.search(r'[0-9]', raw):
            F.add('warn', 'no-number', b.start, min(b.end, b.start + 60),
                  'paragraph in Results and discussion carries no number, atom label or comparison value')


UNIT_RULES = [
    (re.compile(r'\bkcal\s*/\s*mol\b'), 'write kcal mol⁻¹, not kcal/mol'),
    (re.compile(r'\bkJ\s*/\s*mol\b'), 'write kJ mol⁻¹, not kJ/mol'),
    (re.compile(r'\b(?:kcal|kJ)\s*[·.]\s*mol'), 'write "kcal mol⁻¹" (space, no dot)'),
    (re.compile(r'(?<![\w])(cm|mol|eV|Å|nm|mm|pm|g|s)(?:\^\{?\s*[-−–]?\s*\d\s*\}?\^?|[-−–]\d)(?![\w])'),
     'use a Unicode superscript exponent (cm⁻¹, mol⁻¹, eV⁻¹, Å⁻³)'),
    (re.compile(r'Å[23](?![\d])'), 'use Å² / Å³'),
    (re.compile(r'(?<![\w.])[-−+]?\d+\.\d+[    ]?A(?![\w\'′])'), "'A' used for Å"),
    (re.compile(r'\b[Aa]ngstr[oö]ms?\b'), 'use the symbol Å'),
    (re.compile(r'(?<![\w])\d+(?:\.\d+)?\s*(?:deg|degrees?)\b'), 'use ° for angles'),
    (re.compile(r'(?<![\w.])\d+(?:\.\d+)?\s*Debyes?\b'), 'use D for the dipole moment'),
]
GLUED_RE = re.compile(r'(?<![\w.,])([-−+]?\d+(?:\.\d+)?)(Å|eV|kcal|kJ|cm⁻¹|cm|nm|pm|mm|mHa|Ha|a\.u\.|au|Hz|kHz|MHz|GHz|'
                      r'ppm|Da|°C|K|D)(?![\w])')


def check_units(doc, F):
    blocks = body_blocks(doc)
    w = doc.work
    for rx, msg in UNIT_RULES:
        for m in rx.finditer(w):
            if in_blocks(blocks, m.start()):
                F.add('warn', 'units', m.start(), m.end(), msg)
    for m in GLUED_RE.finditer(w):
        if not in_blocks(blocks, m.start()):
            continue
        num, unit = m.group(1), m.group(2)
        if unit == 'D' and '.' not in num:
            continue                 # 2D, 3D
        if unit == 'K' and '.' not in num and len(num.lstrip('-−+')) < 2:
            continue
        F.add('warn', 'units', m.start(), m.end(), f'number glued to its unit: "{m.group()}" (put a space before {unit})')
    au1 = [m for m in re.finditer(r'\d\s*a\.u\.', w) if in_blocks(blocks, m.start())]
    au2 = [m for m in re.finditer(r'\d\s*au\b', w) if in_blocks(blocks, m.start())]
    if au1 and au2:
        F.add('warn', 'units', au2[0].start(), au2[0].end(), 'both "a.u." and "au" used: pick one')


def check_style_extras(doc, F):
    blocks = body_blocks(doc)
    w = doc.work
    for b in blocks:
        if b.back:
            continue
        seg = w[b.start:b.end]
        if b.kind in ('para', 'abstract', 'list', 'quote'):
            for m in re.finditer(r'!(?!\[)', seg):
                F.add('warn', 'style', b.start + m.start(), b.start + m.end(), 'exclamation mark')
            for m in re.finditer(r'\?', seg):
                F.add('warn', 'style', b.start + m.start(), b.start + m.end(), 'question mark (rhetorical question?)')
        for m in re.finditer(r'\*\*[^*\n]+\*\*|(?<!\w)__[^_\n]+__(?!\w)', seg):
            F.add('warn', 'style', b.start + m.start(), b.start + m.end(), 'bold for emphasis')
        if b.kind == 'list':
            F.add('warn', 'style', b.start, b.start + 30, 'bullet/numbered list in the body: write prose')
    hedges = re.compile(r'\b(?:may|might|could|can)\s+(?:potentially|possibly|perhaps|conceivably)\b|'
                        r'\b(?:potentially|possibly|perhaps)\s+(?:suggests?|indicates?|implies?)\b|'
                        r'\bit\s+(?:may|might|could)\s+be\s+(?:possible|suggested)\b|'
                        r'\bseems?\s+to\s+(?:potentially|possibly|suggest)\b', re.I)
    for m in hedges.finditer(w):
        if in_blocks(blocks, m.start()):
            F.add('warn', 'style', m.start(), m.end(), 'stacked hedge')
    for m in re.finditer(r'\bon\s+(?:the\s+)?one\s+hand\b|\bon\s+the\s+other\s+hand\b', w, re.I):
        if in_blocks(blocks, m.start()):
            F.add('warn', 'style', m.start(), m.end(), '"on one hand / on the other hand" framing')
    for m in re.finditer(r'(?<=[\w′])(\.{3}|…|⋯|·{4,}|(?<!·)·{1,2}(?!·))(?=[A-Z])', w):
        if in_blocks(blocks, m.start()):
            F.add('warn', 'style', m.start(), m.end(), 'hydrogen-bond contact: use three middle dots (···)')
    for m in re.finditer(r'\b([A-Z][a-z]?\d*[A-Za-z]?′?)-(H\d*[A-Za-z]?′?|[A-Z][a-z]?\d+[A-Za-z]?′?)(?=[\s·.,;:)]|$)', w):
        if in_blocks(blocks, m.start()) and re.search(r'\d|^[CNOS]$|^H', m.group(1) + ' ' + m.group(2)):
            F.add('warn', 'style', m.start(), m.end(), f'bond "{m.group()}": use an en dash (–)')
    # lengths
    for b in doc.blocks:
        if b.kind == 'abstract':
            n = len(doc.words(b.start, b.end))
            if not 180 <= n <= 250:
                F.add('warn', 'style', b.start, b.start, f'abstract has {n} words (STYLE: 180–250)', ctx=False)
        if b.kind == 'keywords':
            kws = [x for x in doc.work[b.start:b.end].split('\n') if x.strip()]
            if len(kws) == 1:
                kws = [x for x in re.split(r'[;,]', kws[0]) if x.strip()]
            if not 5 <= len(kws) <= 6:
                F.add('warn', 'style', b.start, b.start, f'{len(kws)} keywords (STYLE: 5–6)', ctx=False)
    conc = [b for b in doc.blocks if b.kind in ('para', 'list') and
            strip_number(b.sec_top).lower().startswith('conclusion')]
    if conc:
        n = sum(len(doc.words(b.start, b.end)) for b in conc)
        if not 250 <= n <= 400:
            F.add('warn', 'style', conc[0].start, conc[0].start, f'Conclusions have {n} words (STYLE: 250–400)', ctx=False)


# Docking wording (STYLE.md "Also avoid" and "The compound, its models and its species"). Added in review,
# 28 Sep 2026: parse_style() reads only the banned-word list, so these SAPH rules were written down but not
# checked. Warnings, not HARD: measured K_i values and binding affinities from the literature are legitimate
# (the hCA II rank-order check quotes them); a docking number never is.
DOCK_RULES = [
    (re.compile(r'\bbinding\s+(?:free\s+)?(?:energ(?:y|ies)|affinit(?:y|ies))\b', re.I),
     'binding energy/affinity: only for measured values; a Vina/AutoDock4Zn number is a score (kcal mol⁻¹)'),
    (re.compile(r'\bK~?[id]~?(?![\w~])|\binhibition\s+constants?\b|\bIC~?50~?\b'),
     'K_i/K_d/IC50: never derived from a docking score; quote measured values only, with their source'),
    (re.compile(r'\b(?:strong|stronger|strongest|excellent|good|high|tight|potent)\s+'
                r'(?:binding|binders?|affinit(?:y|ies)|inhibit\w*)\b', re.I),
     'docking cliché ("strong binding", "potent inhibitor"): report the score and the rank against the controls'),
    (re.compile(r'\b(?:drug[- ]candidates?|promising\s+leads?|lead\s+compounds?|potent\s+inhibitors?)\b', re.I),
     'docking cliché: a docking score does not make a drug candidate or a lead'),
    (re.compile(r'\b(?:SAPH|the\s+salt|(?:the\s+)?ion\s+pairs?)\b[^.;]{0,60}?\bdock(?:ed|ing)\b|'
                r'\bdock(?:ed|ing)\s+(?:of\s+)?(?:SAPH|the\s+salt|the\s+ion\s+pairs?)\b', re.I),
     'docked species: never SAPH, the salt or the ion pair; neutral sulfanilamide, the sulfonamidate, the cation'),
]


def check_docking(doc, F):
    blocks = body_blocks(doc)
    w = doc.work
    for rx, msg in DOCK_RULES:
        for m in rx.finditer(w):
            b = in_blocks(blocks, m.start())
            if b is not None and not b.back:
                F.add('warn', 'docking', m.start(), m.end(), msg)


def bigrams(doc, top=30):
    cnt = Counter()
    for b in doc.blocks:
        if b.kind not in ('title', 'abstract', 'para', 'list', 'figure', 'tcaption', 'quote') or b.back or b.refs:
            continue
        seg = doc.work[b.start:b.end].lower()
        for chunk in re.split(r'[.,;:!?()\[\]{}"“”]|\s{3,}', seg):
            toks = re.findall(r"[^\W\d_][\w’'′\-–·]*", chunk)
            for x, y in zip(toks, toks[1:]):
                if x in STOP or y in STOP or len(x) < 2 or len(y) < 2:
                    continue
                cnt[(x, y)] += 1
    return [(f'{x} {y}', n) for (x, y), n in cnt.most_common(top) if n > 1]


# --------------------------------------------------------------------------------------------------------------
ORDER = ['banned', 'em-dash', 'schiff', 'symmetry', 'floats-hard', 'floats', 'openers', 'restatement', 'citations',
         'rhythm', 'no-number', 'units', 'style', 'docking']


def lint(path, rules, top=30):
    doc = Doc.from_file(path)
    F = Findings(doc)
    stats = []
    check_banned(doc, F, rules)
    check_dashes_schiff_symmetry(doc, F)
    check_floats(doc, F)
    check_sentences(doc, F, stats)
    check_citations(doc, F)
    check_numbers_present(doc, F)
    check_units(doc, F)
    check_style_extras(doc, F)
    check_docking(doc, F)
    return doc, F, stats, bigrams(doc, top)


def report(path, doc, F, stats, bgr, quiet=False):
    nwords = sum(len(doc.words(b.start, b.end)) for b in doc.blocks
                 if b.kind in ('abstract', 'para', 'list', 'quote') and not b.back and not b.refs)
    print(f'== lint_prose: {rel(path)}  ({len(doc.line_starts)} lines, {nwords} words of running prose)')
    by = defaultdict(list)
    for it in F.items:
        by[it[1]].append(it)
    for chk in ORDER:
        its = sorted(by.get(chk, []), key=lambda x: x[2])
        if not its or (quiet and all(i[0] != 'HARD' for i in its)):
            continue
        name = {'floats-hard': 'floats (references)', 'floats': 'floats (order/labels)'}.get(chk, chk)
        print(f'\n-- {name} ({len(its)})')
        for sev, _, line, msg, ctx in its:
            if quiet and sev != 'HARD':
                continue
            print(f'  {sev:4s} L{line:<5d} {msg}')
            if ctx:
                print(f'             {ctx}')
    if not quiet:
        print('\n-- sentence length by section (words)')
        print(f'  {"section":46s} {"n":>4s} {"mean":>6s} {"sd":>6s} {"CV":>5s}')
        for sec, n, mean, sd, cv in stats:
            flag = '  LOW' if n >= 5 and cv < 0.35 else ''
            print(f'  {sec[:46]:46s} {n:4d} {mean:6.1f} {sd:6.1f} {cv:5.2f}{flag}')
        print(f'\n-- top content bigrams (count > 1)')
        for i in range(0, len(bgr), 3):
            print('  ' + ''.join(f'{g:34s}{n:3d}   ' for g, n in bgr[i:i + 3]))
    hard = F.hard()
    cnt = Counter(i[1] for i in F.items)
    parts = ', '.join(f'{k.replace("-hard", "")} {v}' for k, v in sorted(cnt.items(), key=lambda kv: ORDER.index(kv[0])))
    print(f'\n-- summary: {len(hard)} hard, {len(F.items) - len(hard)} warnings ({parts or "none"})'
          f' -> {"FAIL" if hard else "PASS"}')
    return len(hard)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('manuscripts', nargs='+', help='manuscript.md path, paper id or paper directory')
    ap.add_argument('--style', default=str(STYLE))
    ap.add_argument('--json', help='also write the findings as JSON')
    ap.add_argument('--quiet', action='store_true', help='print hard violations and the summary only')
    ap.add_argument('--bigrams', type=int, default=30)
    ap.add_argument('--list-rules', action='store_true', help='print the rules parsed from STYLE.md and exit')
    a = ap.parse_args(argv)
    rules = parse_style(a.style)
    if a.list_rules:
        for r in rules:
            print(f'{r.label:40s} limit={r.limit} mode={r.mode:11s} {r.regex.pattern[:110]}')
        return 0
    total_hard, out = 0, {}
    for i, m in enumerate(a.manuscripts):
        p = resolve_manuscript(m)
        if not p.exists():
            print(f'{m}: {p} not found', file=sys.stderr)
            total_hard += 1
            continue
        if i:
            print()
        doc, F, stats, bgr = lint(p, rules, a.bigrams)
        total_hard += report(p, doc, F, stats, bgr, a.quiet)
        out[rel(p)] = {'findings': [dict(zip(('severity', 'check', 'line', 'message', 'context'), it)) for it in F.items],
                       'sentence_stats': [dict(zip(('section', 'n', 'mean', 'sd', 'cv'), s)) for s in stats],
                       'bigrams': bgr}
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding='utf-8')
    return 1 if total_hard else 0


if __name__ == '__main__':
    sys.exit(main())
