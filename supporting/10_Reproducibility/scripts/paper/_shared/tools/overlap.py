#!/usr/bin/env python3
"""overlap.py -- text-similarity check of the manuscripts against each other and against prior text.

Usage
  python3 paper/_shared/tools/overlap.py [--root REPO] [--manuscripts FILE ...] [--allow FILE]
          [--json out.json] [--show N] [--no-abstracts] [--extra-source FILE ...] [--sister DIR | --no-sister]

  By default the manuscripts are paper/main/manuscript.md and paper/si/manuscript.md (those that exist) under
  --root (default: this repository). --manuscripts checks arbitrary files instead (still compared pairwise).
  Standard library only; python-docx is used for .docx text when it is importable (tools/paper-env/bin/python),
  otherwise the same text is read from word/document.xml directly.

Prior text (adapted from cqc-dft to SAPH; paths relative to --root unless stated)
  *.docx in the repo root               the colleague's R4 template ("Revised Manuscript (R4) ...docx")
  PLAN.md, molecules/*/results/*/RESULTS.md, lit/REVIEW.md, lit/notes_*.md (the pattern paper/data/facts.md of the
                                        28 Sep seed matched no file and was removed on 29 Sep 2026)
  lit/prior_text/*.txt|.md|.docx        texts not in the API caches (full texts or abstracts no API returns)
  --sister DIR (default /root/workspace/cqc-dft when it exists): the four published-to-be sister papers
                                        DIR/papers/{3atb,pzc,btpc,btp}/manuscript.md and their si/manuscript.md
  lit/cache/openalex/*.json             title + abstract rebuilt from abstract_inverted_index
  lit/cache/crossref/*.json             abstract (JATS stripped) when OpenAlex has none for that DOI

Method
  Text is lower-cased; citations, cross-references, math, code, comments, link targets, attributes and
  YAML keys are removed (Markdown parsed with mdtext.py); numbers and ranges are dropped, labels such as N4 or
  B3LYP stay, hyphenated words and bond notation (N4–H4···O11) count as one word. Word 7-grams and 10-grams of
  each manuscript are looked up in every other text (n-grams never cross a paragraph, heading or table-cell
  block); runs of overlapping matches are merged into passages and reported with the manuscript line range,
  the source and its line/paragraph. Coverage = % of the manuscript's (checked) words that lie inside a
  7-gram shared with that source. Only windows made entirely of manuscript vocabulary are hashed, so a few
  thousand abstracts take seconds.

Allow-list (paper/_shared/tools/overlap_allow.txt)
  plain line        a boilerplate phrase, normalised like the text; its words are removed from every text
                    before the n-grams are formed (n-grams bridge the gap, so copied text around it still counts)
  re: <regex>       the same, as a regex over the normalised text (lower-case words joined by single spaces)
  !skip-tables      pipe-table rows of the manuscripts are not checked (table captions are)
  !skip-section RX  manuscript sections whose heading (number stripped) matches RX (case-insensitive)
  !exclude-source G prior-text sources whose repo-relative path matches the glob G are ignored
  # ...             comment

Exit status 1 if any 10-gram is shared with a prior-text source, or if a manuscript has more than 3 % of its
words covered by 7-grams shared with another manuscript.
"""
import argparse
import bisect
import fnmatch
import glob
import hashlib
import json
import os
import re
import sys
import time
import zipfile
from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mdtext import Doc, REPO, PAPER_IDS, strip_number  # noqa: E402

HERE = Path(__file__).resolve().parent
DEFAULT_ALLOW = HERE / 'overlap_allow.txt'
N_SHORT, N_LONG = 7, 10
PAIR_LIMIT = 3.0
# a token is a word; hyphenated compounds, bond notation (N4–H4···O11) and elisions stay one token
TOK_RE = re.compile(r"[^\W_]+(?:(?:['’\-‐‑–]|···|⋯|…)[^\W_]+)*")
AU_RE = re.compile(r'\ba\.u\.')
NUM_RE = re.compile(r'[\d.,\-–‐‑]+')
DIGITS = frozenset('0123456789')          # tokens made only of digits and separators (numbers, ranges)


# --------------------------------------------------------------------------------------------------------------
# texts
# --------------------------------------------------------------------------------------------------------------
class Text:
    """A tokenised text. toks: normalised tokens kept for n-grams; segs: paragraph/block id of each token (n-grams
    never cross a boundary); offs/ends: character offsets of each token in `raw` (kept for manuscripts, rebuilt
    on demand for prior text); loc_fn(offset) -> a location label such as 'L12', '¶3' or 'abstract'."""
    __slots__ = ('name', 'kind', 'raw', 'toks', 'offs', 'ends', 'segs', 'prior', 'n_all', 'loc_fn', 'pieces', 'keep')

    def __init__(self, name, kind, prior, loc_fn=None, offsets=False):
        self.name, self.kind, self.prior = name, kind, prior
        self.raw = ''
        self.toks, self.segs = [], []
        self.offs, self.ends = ([], []) if offsets else (None, None)
        self.pieces = []                                 # (text, base, seg) as tokenised
        self.keep = None                                 # token indices kept by the allow-list (lazy offsets)
        self.n_all = 0
        self.loc_fn = loc_fn or (lambda off: '')

    def ensure_offsets(self):
        if self.offs is not None:
            return
        offs, ends = [], []
        for text, base, _ in self.pieces:
            for w, x, y in _tokens(text, True):
                offs.append(base + x)
                ends.append(base + y)
        if self.keep is not None:
            offs = [offs[i] for i in self.keep]
            ends = [ends[i] for i in self.keep]
        self.offs, self.ends = offs, ends

    def loc(self, i):
        self.ensure_offsets()
        return self.loc_fn(self.offs[i])


def _prep(text):
    text = AU_RE.sub('au  ', text)                       # 'a.u.' -> 'au  ' (same length, offsets kept)
    text = text.replace('’', "'").replace('‐', '-').replace('‑', '-')
    low = text.lower()
    return (low, True) if len(low) == len(text) else (text, False)


def _tokens(text, with_offsets):
    """normalised tokens of text (numbers dropped); with_offsets -> [(tok, start, end)]"""
    t, lowered = _prep(text)
    fullnum = NUM_RE.fullmatch
    if with_offsets:
        return [(m.group() if lowered else m.group().lower(), m.start(), m.end())
                for m in TOK_RE.finditer(t) if not (m.group()[0] in DIGITS and fullnum(m.group()))]
    ws = [w for w in TOK_RE.findall(t) if not (w[0] in DIGITS and fullnum(w))]
    return ws if lowered else [w.lower() for w in ws]


def tokenize_into(t, text, base=0, seg=0):
    """append the tokens of `text` (at offset `base` in t.raw, paragraph/block `seg`) to t"""
    t.pieces.append((text, base, seg))
    if t.offs is not None:
        for w, x, y in _tokens(text, True):
            t.toks.append(w)
            t.offs.append(base + x)
            t.ends.append(base + y)
            t.segs.append(seg)
    else:
        ws = _tokens(text, False)
        t.toks.extend(ws)
        t.segs.extend([seg] * len(ws))


class Allow:
    def __init__(self, path):
        self.phrases, self.regexes, self.skip_sections, self.exclude = [], [], [], []
        self.skip_tables = False
        self.path = path
        if path and Path(path).exists():
            for line in Path(path).read_text(encoding='utf-8').splitlines():
                s = line.strip()
                if not s or s.startswith('#'):
                    continue
                if s.startswith('!skip-tables'):
                    self.skip_tables = True
                elif s.startswith('!skip-section'):
                    self.skip_sections.append(re.compile(s.split(None, 1)[1].strip(), re.I))
                elif s.startswith('!exclude-source'):
                    self.exclude.append(s.split(None, 1)[1].strip())
                elif s.startswith('re:'):
                    self.regexes.append(re.compile(s[3:].strip()))
                else:
                    toks = _tokens(s, False)
                    if toks:
                        self.phrases.append(tuple(toks))
        self.first = defaultdict(list)
        for p in self.phrases:
            self.first[p[0]].append(p)

    def apply(self, t):
        """remove allow-listed words from a Text (in place); returns the number removed"""
        n = len(t.toks)
        drop = bytearray(n)
        toks = t.toks
        first = self.first
        for i, w in enumerate(toks):
            if w in first:
                for p in first[w]:
                    L = len(p)
                    if tuple(toks[i:i + L]) == p:
                        for j in range(i, i + L):
                            drop[j] = 1
        if self.regexes:
            starts, pos = [], 0
            for w in toks:
                starts.append(pos)
                pos += len(w) + 1
            s = ' '.join(toks)
            for rx in self.regexes:
                for m in rx.finditer(s):
                    i0 = bisect.bisect_right(starts, m.start()) - 1
                    i1 = bisect.bisect_left(starts, m.end())
                    for j in range(max(0, i0), min(n, i1)):
                        drop[j] = 1
        k = sum(drop)
        if k:
            keep = [i for i in range(n) if not drop[i]]
            for attr in ('toks', 'offs', 'ends', 'segs'):
                v = getattr(t, attr)
                if v is not None:
                    setattr(t, attr, [v[i] for i in keep])
            t.keep = keep if t.keep is None else [t.keep[i] for i in keep]
        return k

    def excluded(self, relpath):
        return any(fnmatch.fnmatch(relpath, g) for g in self.exclude)


def md_text(path, name, kind, prior, allow, manuscript=False):
    doc = Doc.from_file(path)
    t = Text(name, kind, prior, lambda off: f'L{doc.line_of(off)}', offsets=manuscript)
    t.raw = doc.text
    kinds = ('title', 'abstract', 'keywords', 'heading', 'para', 'figure', 'tcaption', 'table', 'list', 'quote')
    for bi, b in enumerate(doc.blocks):
        if b.kind not in kinds or b.refs:
            continue
        if manuscript:
            if allow.skip_tables and b.kind == 'table':
                continue
            if any(rx.search(strip_number(b.sec_top)) or rx.search(strip_number(b.sec_fine))
                   for rx in allow.skip_sections):
                continue
        tokenize_into(t, doc.work[b.start:b.end], b.start, bi)
    t.n_all = len(t.toks)
    return t


def plain_text(raw, name, kind, prior, loc_prefix='L', paragraphs=None, labels=None):
    """paragraphs=None: paragraphs are separated by blank lines and located by line number; otherwise each
    paragraph is a segment located as <loc_prefix><n> (or by labels[n])"""
    if paragraphs is None:
        starts = [0] + [m.end() for m in re.finditer('\n', raw)]
        t = Text(name, kind, prior, lambda off: f'{loc_prefix}{bisect.bisect_right(starts, off)}')
        t.raw = raw
        for k, m in enumerate(re.finditer(r'(?:[^\n]|\n(?![ \t]*\n))+', raw)):
            tokenize_into(t, m.group(), m.start(), k)
    else:
        pstarts, base = [], 0
        for p in paragraphs:
            pstarts.append(base)
            base += len(p) + 1

        def loc(off):
            i = bisect.bisect_right(pstarts, off) - 1
            return labels[i] if labels else f'{loc_prefix}{i + 1}'
        t = Text(name, kind, prior, loc)
        t.raw = '\n'.join(paragraphs)
        for i, p in enumerate(paragraphs):
            tokenize_into(t, p, pstarts[i], i)
    t.n_all = len(t.toks)
    return t


def docx_paragraphs(path):
    try:
        import docx                                          # python-docx (paper-env)
        from docx.text.paragraph import Paragraph
        d = docx.Document(str(path))
        out = [Paragraph(p, d).text for p in d.element.body.iter('{http://schemas.openxmlformats.org/'
                                                                   'wordprocessingml/2006/main}p')]
        return out, 'python-docx'
    except ImportError:
        pass
    W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read('word/document.xml'))
    out = []
    for p in root.iter(W + 'p'):
        buf = []
        for el in p.iter():
            if el.tag == W + 't' and el.text:
                buf.append(el.text)
            elif el.tag == W + 'tab':
                buf.append('\t')
            elif el.tag in (W + 'br', W + 'cr'):
                buf.append('\n')
        out.append(''.join(buf))
    return out, 'zip-xml'


_DEC = json.JSONDecoder()


def _json_field(raw, name):
    """decode only one top-level field of a large JSON record (OpenAlex works are ~40 kB; the abstract is small)"""
    m = re.search(r'"%s"\s*:\s*' % re.escape(name), raw)
    return _DEC.raw_decode(raw, m.end())[0] if m else None


def openalex_text(path):
    raw = Path(path).read_text(encoding='utf-8')
    try:
        d = {k: _json_field(raw, k) for k in ('doi', 'title', 'display_name', 'abstract_inverted_index')}
    except ValueError:
        d = json.loads(raw)
    inv = d.get('abstract_inverted_index') or {}
    words = []
    if inv:
        n = 1 + max(p for ps in inv.values() for p in ps)
        words = [''] * n
        for w, ps in inv.items():
            for p in ps:
                words[p] = w
    doi = (d.get('doi') or '').replace('https://doi.org/', '').lower()
    title = d.get('title') or d.get('display_name') or ''
    return doi, title, ' '.join(w for w in words if w)


def crossref_abstract(path):
    d = json.loads(Path(path).read_text(encoding='utf-8'))
    m = d.get('message', d) if isinstance(d, dict) else {}
    ab = m.get('abstract') if isinstance(m, dict) else None
    if not ab:
        return None, None
    ab = re.sub(r'<[^>]+>', ' ', ab)
    return (m.get('DOI') or '').lower(), re.sub(r'\s+', ' ', ab).strip()


def gather_prior(root, allow, no_abstracts, extra, log, sister=None):
    root = Path(root)
    out = []

    roots = os.path.abspath(root)

    def relp(p):
        r = os.path.relpath(os.path.abspath(p), roots)
        return str(p) if r.startswith('..') else r

    def add(t, path):
        if allow.excluded(relp(path)):
            log['excluded'].append(relp(path))
            return
        allow.apply(t)
        out.append(t)

    # the colleague's R4 template (and any other Word file dropped in the repo root)
    for p in sorted(glob.glob(str(root / '*.docx'))):
        paras, how = docx_paragraphs(p)
        log['docx_reader'] = how
        add(plain_text('\n'.join(paras), relp(p), 'template', True, '¶', paras), p)
    for pat in ('PLAN.md', 'molecules/*/results/*/RESULTS.md', 'lit/REVIEW.md', 'lit/notes_*.md'):
        for p in sorted(glob.glob(str(root / pat))):
            add(md_text(p, relp(p), 'notes', True, allow), p)
    # full texts or abstracts that are not in the API caches: drop them into lit/prior_text/ as .txt, .md or .docx
    for p in sorted(glob.glob(str(root / 'lit/prior_text/*'))):
        if p.endswith('.docx'):
            paras, _ = docx_paragraphs(p)
            add(plain_text('\n'.join(paras), relp(p), 'prior', True, '¶', paras), p)
        elif p.endswith('.md'):
            add(md_text(p, relp(p), 'prior', True, allow), p)
        elif p.endswith('.txt'):
            add(plain_text(Path(p).read_text(encoding='utf-8', errors='replace'), relp(p), 'prior', True), p)
    if sister and Path(sister).is_dir():
        for sid in ('3atb', 'pzc', 'btpc', 'btp'):
            for sub in ('manuscript.md', 'si/manuscript.md'):
                p = Path(sister) / 'papers' / sid / sub
                if p.exists():
                    add(md_text(str(p), f'cqc-dft:{sid}/{sub}', 'sister', True, allow), str(p))
    for p in extra or ():
        if p.endswith('.md'):
            add(md_text(p, relp(p), 'extra', True, allow), p)
        elif p.endswith('.docx'):
            paras, _ = docx_paragraphs(p)
            add(plain_text('\n'.join(paras), relp(p), 'extra', True, '¶', paras), p)
        else:
            add(plain_text(Path(p).read_text(encoding='utf-8', errors='replace'), relp(p), 'extra', True), p)
    if not no_abstracts:
        dois = set()
        for p in sorted(glob.glob(str(root / 'lit/cache/openalex/*.json'))):
            try:
                doi, title, ab = openalex_text(p)
            except (ValueError, OSError):
                log['unreadable'].append(relp(p))
                continue
            if not ab and not title:
                continue
            dois.add(doi)
            t = plain_text('', f'openalex:{doi or Path(p).stem}', 'abstract', True, '', [title, ab],
                           ['title', 'abstract'])
            add(t, p)
            log['abstracts'] += 1
        for p in sorted(glob.glob(str(root / 'lit/cache/crossref/*.json'))):
            try:
                doi, ab = crossref_abstract(p)
            except (ValueError, OSError):
                continue
            if not ab or (doi and doi in dois):
                continue
            add(plain_text('', f'crossref:{doi or Path(p).stem}', 'abstract', True, '', [ab], ['abstract']), p)
            log['abstracts'] += 1
    return out


# --------------------------------------------------------------------------------------------------------------
# matching
# --------------------------------------------------------------------------------------------------------------
def build_index(manuscripts):
    idx = {N_SHORT: defaultdict(list), N_LONG: defaultdict(list)}
    vocab = set()
    for mi, m in enumerate(manuscripts):
        toks, segs = m.toks, m.segs
        vocab.update(toks)
        for n in (N_SHORT, N_LONG):
            d = idx[n]
            for i in range(len(toks) - n + 1):
                if segs[i] == segs[i + n - 1]:
                    d[tuple(toks[i:i + n])].append((mi, i))
    return idx, vocab


def scan(src, idx, vocab, skip_mi=None):
    """-> {mi: {n: [(m_pos, s_pos), ...]}} for every n-gram of src found in the manuscripts"""
    toks, segs = src.toks, src.segs
    L = len(toks)
    res = defaultdict(lambda: {N_SHORT: [], N_LONG: []})
    i = 0
    while i < L:
        if toks[i] not in vocab:
            i += 1
            continue
        j = i
        while j < L and toks[j] in vocab:
            j += 1
        if j - i >= N_SHORT:                         # a run of in-vocabulary words long enough for a 7-gram
            for n in (N_SHORT, N_LONG):
                d = idx[n]
                for k in range(i, j - n + 1):
                    if segs[k] != segs[k + n - 1]:
                        continue
                    hits = d.get(tuple(toks[k:k + n]))
                    if hits:
                        for mi, mp in hits:
                            if mi != skip_mi:
                                res[mi][n].append((mp, k))
        i = j
    return res


def passages(m, src, hits):
    """merge 7-gram hits into passages -> list of dicts, and the set of covered manuscript positions"""
    h7 = sorted(set(hits[N_SHORT]))
    long_starts = sorted({mp for mp, _ in hits[N_LONG]})
    covered = set()
    out = []
    cur = None
    for mp, sp in h7:
        covered.update(range(mp, mp + N_SHORT))
        if cur and mp <= cur['m1']:
            cur['m1'] = max(cur['m1'], mp + N_SHORT)
            continue
        if cur:
            out.append(cur)
        cur = {'m0': mp, 'm1': mp + N_SHORT, 's0': sp}
    if cur:
        out.append(cur)
    for p in out:
        k = bisect.bisect_left(long_starts, p['m0'])
        p['has10'] = k < len(long_starts) and long_starts[k] <= p['m1'] - N_LONG
        p['words'] = p['m1'] - p['m0']
        a = m.offs[p['m0']]
        b = m.ends[p['m1'] - 1]
        p['text'] = re.sub(r'\s+', ' ', m.raw[a:b])
        l0, l1 = m.loc(p['m0']), m.loc(p['m1'] - 1)
        p['where'] = l0 if l0 == l1 else f'{l0}–{l1[1:] if l1.startswith("L") else l1}'
        s0 = p['s0']
        s1 = min(len(src.toks) - 1, s0 + p['words'] - 1)
        sl0, sl1 = src.loc(s0), src.loc(s1)
        p['src_where'] = sl0 if sl0 == sl1 else f'{sl0}–{sl1.lstrip("L¶")}'
    return out, covered


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--root', default=str(REPO), help='repository root (default: this repo)')
    ap.add_argument('--manuscripts', nargs='+', help='check these files instead of paper/{main,si}/manuscript.md')
    ap.add_argument('--sister', default='/root/workspace/cqc-dft',
                    help='cqc-dft checkout whose four manuscripts count as prior text (default: %(default)s)')
    ap.add_argument('--no-sister', action='store_true', help='do not read the sister papers')
    ap.add_argument('--allow', default=None, help='allow-list file (default: overlap_allow.txt next to this script)')
    ap.add_argument('--extra-source', nargs='+', help='additional prior-text files (.md, .docx or plain text)')
    ap.add_argument('--no-abstracts', action='store_true', help='skip the cached OpenAlex/Crossref abstracts')
    ap.add_argument('--show', type=int, default=25, help='passages to print per manuscript/source pair (default 25)')
    ap.add_argument('--json', help='write the full result as JSON')
    a = ap.parse_args(argv)
    t0 = time.time()
    root = Path(a.root)
    allow = Allow(a.allow or DEFAULT_ALLOW)
    if a.manuscripts:
        mpaths = [Path(p) for p in a.manuscripts]
    else:
        mpaths = [root / 'paper' / i / 'manuscript.md' for i in PAPER_IDS if (root / 'paper' / i / 'manuscript.md').exists()]

    def mname(p):
        p = Path(p)
        return p.parent.name if p.name == 'manuscript.md' else p.name
    manus = []
    for p in mpaths:
        t = md_text(p, mname(p), 'manuscript', False, allow, manuscript=True)
        removed = allow.apply(t)
        t.n_all = len(t.toks)
        manus.append((p, t, removed))
    log = defaultdict(list)
    log['abstracts'] = 0
    prior = gather_prior(root, allow, a.no_abstracts, a.extra_source, log, None if a.no_sister else a.sister)
    t1 = time.time()
    print(f'== overlap: {len(manus)} manuscript(s); {len(prior)} prior-text sources '
          f'({log["abstracts"]} abstracts, {sum(len(s.toks) for s in prior):,} words); allow-list '
          f'{allow.path if Path(str(allow.path)).exists() else "(none)"}: {len(allow.phrases)} phrases, '
          f'{len(allow.regexes)} regexes, skip-tables={allow.skip_tables}, {len(allow.skip_sections)} skipped sections')
    if log['duplicates']:
        print('   identical copies read once: ' + '; '.join(log['duplicates']))
    if log.get('docx_reader'):
        print(f'   .docx text via {log["docx_reader"]}')
    if not manus:
        print('   no manuscripts found: nothing to check')
        return 0
    for p, t, removed in manus:
        print(f'   {t.name}: {len(t.toks):,} words checked ({removed} allow-listed words removed)')
    mtexts = [t for _, t, _ in manus]
    idx, vocab = build_index(mtexts)
    fail = []
    result = {'pairs': [], 'prior': []}

    # (a) manuscript pairs
    if len(mtexts) > 1:
        print(f'\n-- (a) manuscript vs manuscript (7-gram coverage of the first by the second; limit {PAIR_LIMIT} %)')
        for sj, src in enumerate(mtexts):
            hits = scan(src, idx, vocab, skip_mi=sj)
            for mi, m in enumerate(mtexts):
                if mi == sj:
                    continue
                h = hits.get(mi, {N_SHORT: [], N_LONG: []})
                ps, cov = passages(m, src, h)
                pct = 100.0 * len(cov) / max(1, len(m.toks))
                bad = pct > PAIR_LIMIT
                n10 = sum(1 for p in ps if p['has10'])
                print(f'  {m.name:>8s} vs {src.name:<8s} {pct:5.2f} %  ({len(cov)} words, {len(ps)} passages, '
                      f'{n10} with a 10-gram){"  FAIL" if bad else ""}')
                if bad:
                    fail.append(f'{m.name} shares {pct:.2f} % with {src.name}')
                for p in sorted(ps, key=lambda p: -p['words'])[:a.show]:
                    print(f'      {"10+" if p["has10"] else " 7+"} {m.name} {p["where"]} ~ {src.name} {p["src_where"]} '
                          f'({p["words"]} w): "{p["text"][:160]}"')
                result['pairs'].append({'manuscript': m.name, 'other': src.name, 'coverage_pct': round(pct, 3),
                                        'words': len(cov), 'passages': ps})

    # (b) prior text
    print('\n-- (b) manuscripts vs prior text (any shared 10-gram fails)')
    per_m = defaultdict(list)
    cov_all = defaultdict(set)
    for src in prior:
        hits = scan(src, idx, vocab)
        for mi, h in hits.items():
            m = mtexts[mi]
            ps, cov = passages(m, src, h)
            if not ps:
                continue
            cov_all[mi] |= cov
            per_m[mi].append((src, ps, cov))
    for mi, m in enumerate(mtexts):
        rows = sorted(per_m.get(mi, []), key=lambda r: -len(r[2]))
        tot = 100.0 * len(cov_all[mi]) / max(1, len(m.toks))
        n10 = sum(1 for _, ps, _ in rows for p in ps if p['has10'])
        print(f'  {m.name}: {tot:.2f} % of words in 7-grams shared with prior text; {len(rows)} sources, '
              f'{n10} passage(s) with a 10-gram{"  FAIL" if n10 else ""}')
        if n10:
            fail.append(f'{m.name}: {n10} passage(s) share a 10-gram with prior text')
        for src, ps, cov in rows:
            pct = 100.0 * len(cov) / max(1, len(m.toks))
            k10 = sum(1 for p in ps if p['has10'])
            print(f'    {pct:5.2f} %  {src.name}  ({len(ps)} passages, {k10} with a 10-gram)')
            for p in sorted(ps, key=lambda p: (-p['has10'], -p['words']))[:a.show]:
                print(f'        {"10+" if p["has10"] else " 7+"} {m.name} {p["where"]} ~ {p["src_where"]} '
                      f'({p["words"]} w): "{p["text"][:160]}"')
            result['prior'].append({'manuscript': m.name, 'source': src.name, 'coverage_pct': round(pct, 3),
                                    'passages': ps})
    print(f'\n-- summary ({time.time() - t0:.1f} s; reading {t1 - t0:.1f} s): '
          + ('FAIL: ' + '; '.join(fail) if fail else 'PASS'))
    if a.json:
        Path(a.json).write_text(json.dumps(result, indent=1, ensure_ascii=False), encoding='utf-8')
    return 1 if fail else 0


if __name__ == '__main__':
    sys.exit(main())
