#!/usr/bin/env python3
"""mdtext.py -- Markdown parsing shared by the manuscript audit tools (standard library only).

Used by lint_prose.py, overlap.py and numbers_audit.py. Not a full Markdown parser: it understands exactly
the subset the manuscripts use (pandoc Markdown, see paper/main/manuscript.md):

  * YAML front matter: title, abstract and keywords are prose; everything else (keys, comments) is masked.
  * Masked (non-prose) spans: HTML comments, fenced code, inline code, $$display$$ and $inline$ math,
    raw HTML tags, link/image targets, {#attribute} blocks, fenced-div lines, citation brackets [@a; @b]
    and bare @keys (cross-references @fig:x, @tbl:x, @eq:x, @sec:x are kept apart from citations).
  * Blocks: title, abstract, keywords, heading, para, figure (captioned image), tcaption ("Table: ..."),
    table (pipe or grid table rows), list, math (a paragraph that is only display math), code.
  * Sections: every block knows its top-level section, its innermost section, whether it is back matter
    (CRediT, declarations, data availability, acknowledgements, ...) and whether it is in the reference list.

Doc.work is the text with every masked character replaced by a space (newlines kept), so offsets and line
numbers in `work` are identical to those in the original text.
"""
import bisect
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent            # paper/_shared/tools
SHARED = HERE.parent                               # paper/_shared
PAPERS = SHARED.parent                             # paper  (main/ = manuscript, si/ = Supporting Information)
REPO = PAPERS.parent
LIT = REPO / 'lit'                                 # the single literature pool (candidates.jsonl) of this project
# Adapted from cqc-dft (four papers, papers/<id>/) to this project's single paper: the ids are the two
# documents, paper/main/manuscript.md and paper/si/manuscript.md.
PAPER_IDS = ('main', 'si')

# mask categories (0 = prose)
PROSE, COMMENT, CODE, MATH, CITE, XREF, ATTR, META, HTML, DIV = range(10)
XREF_KINDS = ('fig', 'tbl', 'eq', 'sec')
PROSE_YAML_KEYS = ('title', 'subtitle', 'abstract', 'keywords')

# back matter (not running prose): Elsevier and ACS headings (ACS Omega: Associated Content, Author Information,
# Acknowledgments, the TOC graphic page "For Table of Contents Only"; added 29 Sep 2026)
BACK_RE = re.compile(r'^(credit|author contributions?|declaration|competing|conflicts? of interest|data availability|'
                     r'acknowledg|funding|generative ai|use of ai|ethic|supplementary|supporting information|'
                     r'associated content|author information|notes$|abbreviations|for table of contents|appendix|'
                     r'references|bibliography|literature cited)', re.I)
REFS_RE = re.compile(r'^(references|bibliography|literature cited)\b', re.I)

KEY_RE = r'[\w][\w:.#$%&\-+?<>~/]*[\w]|[\w]'
CITE_BRACKET_RE = re.compile(r'(?<!!)\[((?:[^\[\]\n]|\n(?![ \t]*\n))*?(?<![\w])-?@(?:[^\[\]\n]|\n(?![ \t]*\n))*?)\]')
CITE_KEY_RE = re.compile(r'(?<![\w@])-?@(' + KEY_RE + ')')
BARE_CITE_RE = re.compile(r'(?<![\w@\\\[])-?@(' + KEY_RE + ')')

# sentence splitting
ABBREV = {'fig', 'figs', 'eq', 'eqs', 'ref', 'refs', 'al', 'vs', 'cf', 'ca', 'approx', 'no', 'nos', 'sect',
          'sects', 'tab', 'resp', 'viz', 'dr', 'prof', 'st', 'vol', 'pp', 'p', 'ed', 'eds', 'e.g', 'i.e',
          'suppl', 'chem', 'phys', 'eqn', 'eqns'}
WORD_RE = re.compile(r"[^\W_]+(?:[’'\-–.·′][^\W_]+)*")


def resolve_manuscript(arg):
    """paper id | directory | file -> Path of the manuscript"""
    if arg in PAPER_IDS:
        return PAPERS / arg / 'manuscript.md'
    p = Path(arg)
    if p.is_dir():
        p = p / 'manuscript.md'
    return p


def rel(p):
    try:
        return str(Path(p).resolve().relative_to(REPO))
    except ValueError:
        return str(p)


def clean_heading(t):
    t = re.sub(r'\{[^{}]*\}\s*$', '', t).strip()
    return t


def strip_number(t):
    return re.sub(r'^\d+(\.\d+)*\.?\s+', '', t).strip()


class Block:
    __slots__ = ('kind', 'start', 'end', 'line0', 'line1', 'sec_top', 'sec_fine', 'back', 'refs', 'level',
                 'title', 'yaml_key')

    def __init__(self, kind, start, end, line0, line1):
        self.kind, self.start, self.end, self.line0, self.line1 = kind, start, end, line0, line1
        self.sec_top = self.sec_fine = ''
        self.back = self.refs = False
        self.level = 0
        self.title = ''
        self.yaml_key = None

    def __repr__(self):
        return f'<{self.kind} L{self.line0}-{self.line1} [{self.sec_fine}]>'


class Citation:
    __slots__ = ('start', 'end', 'keys', 'bracket')

    def __init__(self, start, end, keys, bracket):
        self.start, self.end, self.keys, self.bracket = start, end, keys, bracket

    @property
    def lit_keys(self):
        return [k for k in self.keys if k.split(':', 1)[0] not in XREF_KINDS or ':' not in k]

    @property
    def xref_keys(self):
        return [k for k in self.keys if ':' in k and k.split(':', 1)[0] in XREF_KINDS]


class Doc:
    def __init__(self, text, path=None, markdown=True):
        text = text.replace('\r\n', '\n').replace('\r', '\n')
        self.path = path
        self.text = text
        self.n = len(text)
        self.cat = bytearray(self.n)
        self._work = list(text)
        self.line_starts = [0] + [m.end() for m in re.finditer('\n', text)]
        self.blocks = []
        self.citations = []
        self.eq_labels = []
        self.yaml_end = 0            # offset of the first body character
        if markdown:
            self._parse()
        self.work = ''.join(self._work)

    # -- helpers ------------------------------------------------------------------------------------------
    @classmethod
    def from_file(cls, path, markdown=True):
        return cls(Path(path).read_text(encoding='utf-8', errors='replace'), path=path, markdown=markdown)

    def line_of(self, off):
        return bisect.bisect_right(self.line_starts, off)          # 1-based

    def line_text(self, lineno):
        a = self.line_starts[lineno - 1]
        b = self.line_starts[lineno] - 1 if lineno < len(self.line_starts) else self.n
        return self.text[a:b]

    def context(self, a, b, width=40):
        s = max(0, a - width)
        e = min(self.n, b + width)
        pre = self.text[s:a].replace('\n', ' ')
        mid = self.text[a:b].replace('\n', ' ')
        post = self.text[b:e].replace('\n', ' ')
        return ('…' if s > 0 else '') + pre + '«' + mid + '»' + post + ('…' if e < self.n else '')

    def mark(self, a, b, code):
        cat, work, text = self.cat, self._work, self.text
        for i in range(max(0, a), min(b, self.n)):
            if cat[i] == 0:
                cat[i] = code
                work[i] = '\n' if text[i] == '\n' else ' '

    def _cur(self):
        return ''.join(self._work)

    def masked(self, a, b):
        return self.work[a:b]

    def block_at(self, off):
        for b in self.blocks:
            if b.start <= off < b.end:
                return b
        return None

    # -- parsing ------------------------------------------------------------------------------------------
    def _parse(self):
        text = self.text
        lines = text.split('\n')
        starts = self.line_starts
        # HTML comments
        for m in re.finditer(r'<!--.*?-->', text, re.S):
            self.mark(m.start(), m.end(), COMMENT)
        body_line = self._yaml(lines)
        # fenced code blocks (line based)
        in_code, fence = False, ''
        code_lines = set()
        for i in range(body_line, len(lines)):
            s = lines[i]
            m = re.match(r'^\s*(`{3,}|~{3,})', s)
            if in_code:
                code_lines.add(i)
                if m and m.group(1)[0] == fence[0] and len(m.group(1)) >= len(fence):
                    in_code = False
                continue
            if m:
                in_code, fence = True, m.group(1)
                code_lines.add(i)
        for i in code_lines:
            a = starts[i]
            self.mark(a, a + len(lines[i]) + 1, CODE)
        # fenced divs
        for i in range(body_line, len(lines)):
            if i not in code_lines and re.match(r'^\s*:{3,}', lines[i]):
                a = starts[i]
                self.mark(a, a + len(lines[i]), DIV)
        cur = self._cur()
        # display math (+ equation label)
        for m in re.finditer(r'\$\$(.+?)\$\$([ \t]*\{#(eq:[^}\s]*)[^}]*\})?', cur, re.S):
            self.mark(m.start(), m.end(), MATH)
            if m.group(3):
                self.eq_labels.append((m.group(3), m.start()))
        cur = self._cur()
        for m in re.finditer(r'(`+)(.+?)\1', cur, re.S):
            if '\n\n' not in m.group(2):
                self.mark(m.start(), m.end(), CODE)
        cur = self._cur()
        for m in re.finditer(r'(?<![\\$])\$(?![\s$])((?:[^$\\\n]|\\.|\n(?![ \t]*\n))+?)(?<![\s\\])\$(?!\d)', cur):
            self.mark(m.start(), m.end(), MATH)
        cur = self._cur()
        for m in re.finditer(r'</?[A-Za-z][A-Za-z0-9]*(?:\s[^<>\n]*)?/?>', cur):
            self.mark(m.start(), m.end(), HTML)
        cur = self._cur()
        # link / image targets and attribute blocks
        for m in re.finditer(r'\]\((?:[^()\s]|\([^()\s]*\))*(?:\s+"[^"]*")?\)', cur):
            self.mark(m.start() + 1, m.end(), ATTR)
        for m in re.finditer(r'\{[^{}\n]*\}', cur):
            inner = m.group()[1:-1].strip()
            if inner.startswith(('#', '.', '-')) or '=' in inner:
                self.mark(m.start(), m.end(), ATTR)
        cur = self._cur()
        # citations and cross-references
        for m in CITE_BRACKET_RE.finditer(cur):
            keys = [k.group(1) for k in CITE_KEY_RE.finditer(m.group(1))]
            if not keys:
                continue
            c = Citation(m.start(), m.end(), keys, True)
            self.citations.append(c)
            self.mark(m.start(), m.end(), XREF if not c.lit_keys else CITE)
        cur = self._cur()
        for m in BARE_CITE_RE.finditer(cur):
            if m.start() > 0 and cur[m.start() - 1] == '@':
                continue
            c = Citation(m.start(), m.end(), [m.group(1)], False)
            self.citations.append(c)
            self.mark(m.start(), m.end(), XREF if not c.lit_keys else CITE)
        self.citations.sort(key=lambda c: c.start)
        self.work = self._cur()
        self._blocks(lines, body_line, code_lines)

    def _yaml(self, lines):
        """Mask YAML metadata; create title/abstract/keywords blocks. -> first body line index"""
        if not lines or lines[0].strip() != '---':
            return 0
        end = None
        for i in range(1, len(lines)):
            if lines[i].strip() in ('---', '...'):
                end = i
                break
        if end is None:
            return 0
        starts = self.line_starts
        for i in (0, end):
            self.mark(starts[i], starts[i] + len(lines[i]), META)
        key = None
        acc = {}                        # key -> [first_off, last_off, line0, line1]
        for i in range(1, end):
            s = lines[i]
            a = starts[i]
            if not s.strip() or re.match(r'^\s*#', s):
                self.mark(a, a + len(s), META)
                continue
            m = re.match(r'^([A-Za-z_][\w-]*):(\s*)(.*)$', s)
            if m:
                key = m.group(1).lower()
                val = m.group(3)
                if key in PROSE_YAML_KEYS and val.strip() not in ('', '|', '>', '|-', '>-', '|+', '>+'):
                    self.mark(a, a + m.start(3), META)
                    v0, v1 = a + m.start(3), a + len(s)
                    if val[:1] in '"\'':
                        self.mark(v0, v0 + 1, META)
                    if val.rstrip()[-1:] in '"\'' and len(val.strip()) > 1:
                        e = a + len(s.rstrip()) - 1
                        self.mark(e, e + 1, META)
                    acc.setdefault(key, [v0, v1, i + 1, i + 1])
                    acc[key][1], acc[key][3] = v1, i + 1
                else:
                    self.mark(a, a + len(s), META)
                continue
            if key in PROSE_YAML_KEYS:
                m2 = re.match(r'^(\s*-\s+)(.*)$', s)
                if m2:
                    self.mark(a, a + len(m2.group(1)), META)
                v0 = a + (len(m2.group(1)) if m2 else len(s) - len(s.lstrip()))
                acc.setdefault(key, [v0, a + len(s), i + 1, i + 1])
                acc[key][1], acc[key][3] = a + len(s), i + 1
            else:
                self.mark(a, a + len(s), META)
        for k, (v0, v1, l0, l1) in acc.items():
            b = Block({'subtitle': 'title'}.get(k, k), v0, v1, l0, l1)
            b.yaml_key = k
            b.sec_top = b.sec_fine = {'title': 'Title', 'subtitle': 'Title', 'abstract': 'Abstract',
                                      'keywords': 'Keywords'}[k]
            self.blocks.append(b)
        self.yaml_end = starts[end + 1] if end + 1 < len(starts) else self.n
        return end + 1

    def _blocks(self, lines, first, code_lines):
        starts = self.line_starts
        work_lines = self.work.split('\n')
        sec_top, sec_fine, back, refs = '', '', False, False
        sec_levels = {}
        cur = None                       # [kind, first_line_idx, last_line_idx]
        out = []

        def flush():
            nonlocal cur
            if cur is None:
                return
            kind, i0, i1 = cur
            a = starts[i0]
            b = starts[i1] + len(lines[i1])
            if kind == 'para':
                body = self.work[a:b]
                first_line = lines[i0].lstrip()
                if not body.strip():
                    kind = 'math' if '$$' in self.text[a:b] else 'meta'
                elif first_line.startswith('!['):
                    kind = 'figure'
                elif re.match(r'^([-*+]|\d+[.)])\s+', first_line):
                    kind = 'list'
                elif first_line.startswith('>'):
                    kind = 'quote'
            out.append((kind, a, b, i0 + 1, i1 + 1))
            cur = None

        for i in range(first, len(lines)):
            s = lines[i]
            if i in code_lines:
                if cur and cur[0] != 'code':
                    flush()
                if cur is None:
                    cur = ['code', i, i]
                else:
                    cur[2] = i
                continue
            if not s.strip() or re.match(r'^\s*:{3,}', s):
                flush()
                if re.match(r'^\s*:{3,}\s*\{#refs', s):
                    out.append(('refsdiv', starts[i], starts[i] + len(s), i + 1, i + 1))
                continue
            hm = re.match(r'^(#{1,6})\s+(.*?)\s*$', s)
            if hm:
                flush()
                out.append(('heading', starts[i], starts[i] + len(s), i + 1, i + 1))
                continue
            st = s.lstrip()
            is_table = st.startswith('|') or bool(re.match(r'^\+[-=+:]+\+\s*$', st))
            is_tcap = bool(re.match(r'^(Table:|:)\s', st))
            kind = 'table' if is_table else ('tcaption' if is_tcap else 'para')
            if cur is None:
                cur = [kind, i, i]
            elif (cur[0] == 'table') != is_table or is_tcap or cur[0] == 'code':
                flush()
                cur = [kind, i, i]
            else:
                cur[2] = i
        flush()
        for kind, a, b, l0, l1 in out:
            blk = Block(kind, a, b, l0, l1)
            if kind == 'heading':
                m = re.match(r'^(#{1,6})\s+(.*?)\s*$', self.text[a:b])
                level = len(m.group(1))
                title = clean_heading(m.group(2))
                blk.level, blk.title = level, title
                bare = strip_number(title)
                sec_levels = {k: v for k, v in sec_levels.items() if k < level}
                sec_levels[level] = title
                if level == 1 or not sec_top:
                    sec_top = title
                    back = bool(BACK_RE.match(bare))
                    if REFS_RE.match(bare):
                        refs = True
                sec_fine = title
            if kind == 'refsdiv':
                refs = True
            blk.sec_top, blk.sec_fine, blk.back, blk.refs = sec_top, sec_fine, back, refs
            self.blocks.append(blk)
        self.blocks.sort(key=lambda b: b.start)

    # -- derived views ------------------------------------------------------------------------------------
    def prose_blocks(self, kinds=('title', 'abstract', 'keywords', 'heading', 'para', 'figure', 'tcaption',
                                  'table', 'list', 'quote'), refs=False):
        return [b for b in self.blocks if b.kind in kinds and (refs or not b.refs)]

    def sentences(self, blk):
        """-> list of (start, end) offsets of the sentences in a block (on the masked text)."""
        a, b = blk.start, blk.end
        s = self.work[a:b]
        text = self.text
        out = []
        start = 0
        n = len(s)
        for m in re.finditer(r'[.!?]+', s):
            e = m.end()
            if e < n and s[e] not in ' \t\n   "”’)]*_':
                continue
            j = e
            while j < n and s[j] in '"”’)]*_':
                j += 1
            k = j
            while k < n and s[k] in ' \t\n   ':
                k += 1
            if k >= n:
                break
            if k == j and j == e:
                continue
            nxt = s[k]
            # first visible character of the original text after the punctuation (math is masked in s)
            k0 = e
            while k0 < n and text[a + k0] in ' \t\n   "”’)]*_':
                k0 += 1
            orig_next = text[a + k0] if k0 < n else ''
            if not (nxt.isupper() or nxt in '“"(' or orig_next == '$' or nxt in '*'):
                continue
            if m.group() == '.':
                before = s[start:m.start()].rstrip()
                tok = re.search(r'([A-Za-z]+(?:\.[A-Za-z]+)*)$', before)
                if tok and before.endswith(tok.group(1)):
                    w = tok.group(1)
                    if w.lower() in ABBREV or (len(w) == 1 and w.isupper()):
                        continue
            out.append((a + start, a + e))
            start = k
        if s[start:].strip():
            out.append((a + start, a + len(s.rstrip())))
        # trim leading whitespace
        res = []
        for x, y in out:
            while x < y and self.work[x] in ' \t\n':
                x += 1
            if x < y:
                res.append((x, y))
        return res

    def words(self, a, b):
        return WORD_RE.findall(self.work[a:b])

    # -- figures and tables -------------------------------------------------------------------------------
    def float_defs(self):
        """-> {'fig': [(label|None, block)], 'tbl': [...]} in document order"""
        items = []
        for b in self.blocks:
            raw = self.text[b.start:b.end]
            if b.kind == 'figure':
                cap = re.match(r'!\[(.*?)\]\(', raw, re.S)
                if cap and cap.group(1).strip():
                    lab = re.search(r'\{#(fig:[^}\s]+)', raw)
                    items.append((b.start, 'fig', lab.group(1) if lab else None, b))
            elif b.kind == 'tcaption':
                lab = re.search(r'\{#(tbl:[^}\s]+)', raw)
                items.append((b.start, 'tbl', lab.group(1) if lab else None, b))
        # <!-- include: tables/x.md --> lines (include.lua) define the floats of the included file, in order,
        # at the position of the include line. Paths are relative to the manuscript's directory.
        if self.path is not None:
            base = Path(self.path).resolve().parent
            for m in re.finditer(r'^<!--\s*include:\s*(.+?)\s*-->[ \t]*$', self.text, re.M):
                f = base / m.group(1)
                if not f.is_file():
                    continue
                inc = f.read_text(encoding='utf-8', errors='replace')
                ln = self.line_of(m.start())
                blk = Block('include', m.start(), m.end(), ln, ln)
                for lm in re.finditer(r'\{#((fig|tbl):[^}\s]+)', inc):
                    items.append((m.start(), lm.group(2), lm.group(1), blk))
        items.sort(key=lambda t: t[0])
        figs = [(lab, b) for _, k, lab, b in items if k == 'fig']
        tbls = [(lab, b) for _, k, lab, b in items if k == 'tbl']
        return {'fig': figs, 'tbl': tbls}

    LIT_FLOAT_RE = re.compile(r'\b(Figs?\.|Figures?|Tables?)(?:[   ]+|(?=\d))'
                              r'(S?\d+[a-z]?(?:\s*(?:,|and|&|–|-|to)\s*S?\d+[a-z]?)*)')

    def float_refs(self):
        """-> list of (kind, number:int|None, label|None, offset, literal:bool). Literal 'Fig. 2' text and
        pandoc cross-references (@fig:x) both count."""
        out = []
        for c in self.citations:
            for k in c.xref_keys:
                kind = k.split(':', 1)[0]
                if kind in ('fig', 'tbl'):
                    out.append((kind, None, k, c.start, False))
        # A Supporting Information file (YAML number-prefix: S) numbers its own floats S1, S2, ...; there a
        # literal 'Fig. S3' is one of its own figures and a plain 'Fig. 3' belongs to the main text. In a main
        # text it is the other way round.
        own_s = bool(re.search(r'^number-prefix:\s*["\']?S', self.text[:max(self.yaml_end, 0)], re.M))
        for m in self.LIT_FLOAT_RE.finditer(self.work):
            kind = 'fig' if m.group(1).lower().startswith('fig') else 'tbl'
            nums = []
            parts = re.findall(r'(S?)(\d+)[a-z]?|(–|-|to)', m.group(2))
            prev = None
            rng = False
            for sup, d, dash in parts:
                if dash:
                    rng = True
                    continue
                if bool(sup) != own_s:
                    prev = None
                    continue
                v = int(d)
                if rng and prev is not None and v > prev:
                    nums.extend(range(prev + 1, v + 1))
                else:
                    nums.append(v)
                prev, rng = v, False
            for v in nums:
                out.append((kind, v, None, m.start(), True))
        out.sort(key=lambda r: r[3])
        return out


def ctx_line(doc, a, b, width=40):
    return f'L{doc.line_of(a)}', doc.context(a, b, width)
