#!/usr/bin/env python3
"""Build the SAPH manuscript or SI: paper/<main|si>/manuscript.md -> SAPH_manuscript.docx/.pdf or
SAPH_Supplementary_Material.docx/.pdf (SAPH_Supporting_Information with the ACS profile) in the same directory.

  tools/paper-env/bin/python paper/_shared/tools/build.py main|si|<dir> [--docx] [--pdf] [--spacing 1.5|2]
                                                           [--journal jics|elsevier|acs-omega]
  main = paper/main, si = paper/si, or any directory holding manuscript.md (its name gives the output stem).
  --docx / --pdf select one output (default: both).  --spacing 2 gives the double-spaced submission version.
  --no-bib skips the lit.py step (offline rebuild from an existing references.json).
  --journal selects the citation profile (JOURNALS below). Default jics (Journal of the Indian Chemical Society, the
  authors' decision of 5 Oct 2026): the Elsevier numeric style with titles, "[1-3]" in the text and "[1] A. Author,
  Title, Abbrev. J. Vol (Year) pages. https://doi.org/..." in the list, "Fig. 1" and "Eq. (1)". acs-omega (target
  until 5 Oct 2026): American Chemical Society CSL, superscript citation numbers collapsed to ranges ("salt.^1-3^",
  placed after a following full stop or comma as ACS prints them) and "(1) Author, A.; ... Abbrev. J. Year, Vol
  (Issue), pages." entries.
  (Plain python3 works too: without lxml the script re-executes itself with tools/paper-env/bin/python.)

Adapted from cqc-dft (papers/_shared/tools/build.py, four papers) to this project's single paper; the
toolchain is found at repo-relative paths under the git-ignored tools/ (on the Pi these are symlinks into
cqc-dft/tools, made by paper/_shared/tools/link_toolchain.sh): tools/bin/pandoc, tools/paper-env,
tools/texlive-extracted (luaotfload) and tools/fonts/extracted (TeX Gyre Termes, Liberation).

Steps
  1. python3 paper/_shared/tools/lit.py bib --manuscript <dir>/manuscript.md --out <dir>/references.json
     (CSL-JSON from Crossref, only the keys cited). Unknown citation keys stop the build.
  2. references.json -> build/references.csl.json: HTML entities and tags cleaned, ISO 4 journal
     abbreviations from paper/_shared/csl/journal_abbreviations.json, volume year (print issue) instead of the
     online-first year, article numbers used as pages, Crossref's type "journal-article" mapped to the CSL type
     "article-journal" (the ACS style formats only that type as a journal article). Journals without an
     abbreviation are listed.
  3. DOCX: pandoc + crossref.lua + citeproc (the profile's CSL) + bibfix.lua with paper/_shared/reference.docx, then
     post-processing of word/document.xml: continuous line numbers, three-rule tables with 10 pt "Table Text",
     borderless equation tables; docProps/custom.xml emptied (pandoc writes the build metadata there, local paths
     included). Checked with docx_check.py.
  4. PDF: pandoc --pdf-engine=lualatex with manuscript.latex (TeX Gyre Termes by file path, unicode-math, lineno,
     booktabs). Any "Missing character" stops the build. The .tex and the full LaTeX log stay in build/.
Everything intermediate goes to <dir>/build/ (git-ignored).
"""
import argparse
import html
import json
import os
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent            # paper/_shared/tools
SHARED = HERE.parent                               # paper/_shared
PAPERS = SHARED.parent                             # paper
REPO = PAPERS.parent
PAPER_IDS = ('main', 'si')
STEM = 'SAPH'                                      # output file stem: SAPH_manuscript, SAPH_Supporting_Information
CSL_DIR = SHARED / 'csl'
# Citation profiles. 'citations' tells docx_check.py how the numbers print: 'superscript' (ACS; the build wraps
# each citation in the character style "Citation Number", see bibfix.lua) or 'brackets' ("[1-3]", Elsevier).
# 'xref' is the cross-reference style of crossref.lua: 'acs' prints "Figure 1", "Figure 1." and "eq 1" (the ACS
# journals' labels), 'default' prints "Fig. 1", "Fig. 1." and "Eq. (1)".
# 'si' is the file-name stem of the second document: Elsevier calls it Supplementary material, ACS Supporting
# Information.
JOURNALS = {
    'jics': {'csl': CSL_DIR / 'elsevier-with-titles.csl', 'citations': 'brackets', 'xref': 'default', 'abbrev': 'iso4',
             'si': '_Supplementary_Material',
             'name': 'Journal of the Indian Chemical Society (Elsevier numeric style with titles)'},
    'elsevier': {'csl': CSL_DIR / 'elsevier-with-titles.csl', 'citations': 'brackets', 'xref': 'default', 'abbrev': 'iso4',
                 'si': '_Supplementary_Material', 'name': 'Elsevier numeric style with titles'},
    'acs-omega': {'csl': CSL_DIR / 'american-chemical-society.csl', 'citations': 'superscript', 'xref': 'acs', 'abbrev': 'cassi',
                  'si': '_Supporting_Information', 'name': 'ACS Omega (American Chemical Society style)'},
}
DEFAULT_JOURNAL = 'jics'                           # target journal since 5 Oct 2026 (authors' decision)
SI_STEM = JOURNALS[DEFAULT_JOURNAL]['si']           # set from --journal in main()
CSL = JOURNALS[DEFAULT_JOURNAL]['csl']              # set from --journal in main()
XREF = JOURNALS[DEFAULT_JOURNAL]['xref']            # set from --journal in main()
ABBREV_STYLE = JOURNALS[DEFAULT_JOURNAL]['abbrev']  # set from --journal in main()
ABBREV = SHARED / 'csl' / 'journal_abbreviations.json'
REFDOC = SHARED / 'reference.docx'
TEMPLATE = HERE / 'manuscript.latex'
FILTER = HERE / 'crossref.lua'
INCLUDE = HERE / 'include.lua'
BIBFIX = HERE / 'bibfix.lua'
LIT = HERE / 'lit.py'
PAPER_ENV_PY = REPO / 'tools' / 'paper-env' / 'bin' / 'python'


# -------------------------------------------------------------------------------------------------
# small helpers
# -------------------------------------------------------------------------------------------------
class BuildError(Exception):
    pass


def say(msg):
    print(msg, flush=True)


def rel(p):
    """Path relative to the repository when it lies inside it (for messages)."""
    try:
        return Path(p).resolve().relative_to(REPO)
    except ValueError:
        return Path(p)


def pandoc_bin():
    p = REPO / 'tools' / 'bin' / 'pandoc'
    if p.exists():
        return str(p)
    w = shutil.which('pandoc')
    if not w:
        raise BuildError('pandoc not found (expected tools/bin/pandoc)')
    return w


def run(cmd, cwd=None, env=None, check=True):
    r = subprocess.run([str(c) for c in cmd], cwd=cwd, env=env, capture_output=True, text=True)
    if check and r.returncode != 0:
        raise BuildError(f'command failed ({r.returncode}): {" ".join(map(str, cmd))}\n{r.stdout}\n{r.stderr}')
    return r


def out_stem(ID):
    """SAPH -> SAPH_manuscript ; SAPH_SI -> SAPH_Supplementary_Material (Elsevier) or SAPH_Supporting_Information
    (ACS), by the journal profile"""
    return ID[:-3] + SI_STEM if ID.endswith('_SI') else ID + '_manuscript'


def resolve(paper):
    """-> (paper_dir, ID). Every document of this project cites from the one pool, lit/candidates.jsonl."""
    d = PAPERS / paper if paper in PAPER_IDS else Path(paper).resolve()
    if (d / 'manuscript.md').exists():
        if d.resolve() == (PAPERS / 'main').resolve():
            return d, STEM
        if d.resolve() == (PAPERS / 'si').resolve():
            return d, STEM + '_SI'
        return d, d.name.upper().replace('-', '_')
    raise BuildError(f'unknown paper {paper!r}: use main, si, or a directory holding manuscript.md')


# -------------------------------------------------------------------------------------------------
# 1. references
# -------------------------------------------------------------------------------------------------
def bibliography(paper_dir):
    cmd = ['python3', LIT, 'bib', '--manuscript', str(paper_dir / 'manuscript.md'),
           '--out', str(paper_dir / 'references.json')]
    r = run(cmd, cwd=REPO, check=False)
    out = (r.stdout + r.stderr).strip()
    say('  lit.py bib: ' + out.replace('\n', '\n  '))
    if r.returncode != 0:
        raise BuildError('lit.py bib failed')
    if 'not a verified candidate' in out:
        raise BuildError('the manuscript cites keys that lit.py does not know (see above); add them with '
                         '"lit.py add" or fix the key')


def _norm_title(t):
    t = html.unescape(t or '')
    t = t.replace('–', '-').replace('—', '-').replace('&', 'and')
    t = re.sub(r'[^a-z0-9]+', ' ', t.lower()).strip()
    return re.sub(r'^the ', '', t)


ALLOWED_TAGS = re.compile(r'</?(i|b|sup|sub)>', re.I)


def _clean_markup(s):
    """Crossref titles: keep <i>, <b>, <sup>, <sub> (citeproc understands them), drop everything else."""
    if not isinstance(s, str):
        return s
    s = html.unescape(s)
    s = re.sub(r'<scp>(.*?)</scp>', r'\1', s, flags=re.I | re.S)
    s = re.sub(r'<(?!/?(?:i|b|sup|sub)>|span class="nocase">|/span>)[^>]+>', '', s)   # citeproc reads the nocase span
    s = re.sub('[\u2010\u2011]', '-', s)                    # Wiley deposits U+2010 in titles ('Hydrogen‐Bond')
    return re.sub(r'\s+', ' ', s).strip()


def _year_of(d):
    """Citation year: the print issue (volume year) when there is one, else the earliest publication date."""
    ji = d.get('journal-issue') or {}
    for v in (ji.get('published-print'), d.get('published-print'), d.get('issued'), d.get('published-online'),
              d.get('published')):
        if v and v.get('date-parts') and v['date-parts'][0] and v['date-parts'][0][0]:
            return v['date-parts'][0]
    return None


# Journal abbreviations. The table (journal_abbreviations.json) holds the CASSI forms that ACS journals print. Elsevier
# asks for the List of Title Word Abbreviations (ISO 4): no commas, colons or dashes between the words, and a few
# words abbreviated differently. Irregular cases are listed under "_iso4" in the table.
ISO4_WORDS = {'Delivery': 'Deliv.', 'Discovery': 'Discov.', 'Graphics': 'Graph.', 'Modell.': 'Model.',
              'Cheminf.': 'Cheminform.', 'Kristallogr.': 'Krist.'}


def iso4(ab, special):
    """ISO 4 (LTWA) form of a CASSI abbreviation of the table."""
    if ab in special:
        return special[ab]
    ab = re.sub(r'\s+[-–]\s+', ' ', ab)                     # Chem. - Eur. J.
    ab = re.sub(r'[,:]', '', ab)                            # Acta Crystallogr., Sect. B: Struct. Sci.
    return ' '.join(ISO4_WORDS.get(w, w) for w in ab.split())


KEEP = {'id', 'type', 'title', 'container-title', 'container-title-short', 'author', 'editor', 'issued',
        'volume', 'issue', 'page', 'DOI', 'URL', 'publisher', 'publisher-place', 'edition', 'number',
        'collection-title', 'event-title'}


def normalise_refs(paper_dir, build_dir):
    src = paper_dir / 'references.json'
    if not src.exists():
        raise BuildError(f'{src} missing (run lit.py bib)')
    refs = json.loads(src.read_text())
    abbrev = json.loads(ABBREV.read_text()) if ABBREV.exists() else {}
    amap = {_norm_title(k): v for k, v in abbrev.items() if not k.startswith('_')}
    special = abbrev.get('_iso4', {})
    missing = {}
    out = []
    for d in refs:
        e = {k: v for k, v in d.items() if k in KEEP}
        # Crossref's "journal-article" is not a CSL type; the ACS style formats only "article-journal" as a journal
        # article (abbreviated title, bold year, italic volume, issue), everything else falls to its generic branch.
        # The Elsevier style has no journal-article branch, so its output is unchanged (ACS adaptation, 29 Sep 2026).
        if e.get('type') == 'journal-article':
            e['type'] = 'article-journal'
        for k in ('title', 'container-title', 'collection-title', 'event-title'):
            v = e.get(k)
            if isinstance(v, list):
                v = v[0] if v else ''
            if v:
                e[k] = _clean_markup(v)
        ct = html.unescape(d.get('container-title') if isinstance(d.get('container-title'), str)
                           else ((d.get('container-title') or [''])[0]))
        ct = re.sub(r'\s+', ' ', ct).strip()
        if ct:
            ab = amap.get(_norm_title(ct))
            if ab is None:
                missing[ct] = d.get('container-title-short')
                ab = ct                                  # full title rather than Crossref's mixed-style short
            elif ABBREV_STYLE == 'iso4':
                ab = iso4(ab, special)
            e['container-title-short'] = _clean_markup(ab)
        # Wiley deposits hyphenated given names with U+2010 ('Jan‐Michael'), which citeproc does not split, so
        # 'J.-M. Mewes' printed as 'J. Mewes': use an ASCII hyphen in given names (refs audit, 29 Sep 2026)
        for role in ('author', 'editor'):
            if isinstance(e.get(role), list):
                e[role] = [dict(a, **{k: re.sub('[‐‑]', '-', a[k]) for k in ('given', 'family')
                                      if isinstance(a.get(k), str)}) for a in e[role]]
        y = _year_of(d)
        if y:
            e['issued'] = {'date-parts': [y]}
        if not e.get('page'):
            an = d.get('article-number') or d.get('number')
            if an:
                e['page'] = str(an)
        if e.get('page'):
            e['page'] = e['page'].replace('--', '-')
            # one-page articles that Crossref stores as a range (Acta Cryst. E 'o1488-o1488'): print the page once
            # (refs audit, 29 Sep 2026)
            a, sep, b = e['page'].partition('-')
            if sep and a.strip() == b.strip():
                e['page'] = a.strip()
        out.append(e)
    dst = build_dir / 'references.csl.json'
    dst.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    if missing:
        say(f'  WARNING: {len(missing)} journal(s) without an abbreviation in '
            f'{rel(ABBREV)}; the full title is printed:')
        for k, v in sorted(missing.items()):
            say(f'    "{k}"   (Crossref short form: {v!r})')
    return dst


# -------------------------------------------------------------------------------------------------
# 2. pandoc runs
# -------------------------------------------------------------------------------------------------
def pandoc_common(paper_dir, build_dir, refs):
    return [pandoc_bin(), 'manuscript.md', '--from', 'markdown',
            '--resource-path', str(paper_dir),
            '--lua-filter', str(INCLUDE),
            '--lua-filter', str(FILTER),
            '--citeproc', '--csl', str(CSL), '--bibliography', str(refs),
            '--lua-filter', str(BIBFIX),
            '--metadata', f'crossref-map={build_dir / "crossref.json"}',
            '--metadata', 'link-bibliography=true',
            '--metadata', f'xref-style={XREF}']


FATAL_PATTERNS = [
    (re.compile(r'Citeproc: citation (\S+) not found'), 'citation key not in references.json'),
    (re.compile(r'Could not fetch resource'), 'missing image or resource'),
    (re.compile(r'\[crossref\] ERROR'), 'cross-reference error'),
    (re.compile(r'Missing character'), 'character missing from every font'),
    (re.compile(r'Duplicate (note )?reference', re.I), 'duplicate reference'),
]


def check_messages(stderr, what):
    problems = []
    for line in stderr.splitlines():
        for pat, why in FATAL_PATTERNS:
            if pat.search(line):
                problems.append(f'{why}: {line.strip()}')
    other = [l for l in stderr.splitlines() if l.strip() and not any(p.search(l) for p, _ in FATAL_PATTERNS)]
    for l in other:
        say(f'  [{what}] {l}')
    if problems:
        raise BuildError(f'{what}: ' + '\n  '.join(problems))


def reference_doc(build_dir, spacing):
    if abs(spacing - 1.5) < 1e-9 and REFDOC.exists():
        return REFDOC
    sys.path.insert(0, str(HERE))
    import make_reference_docx
    target = REFDOC if abs(spacing - 1.5) < 1e-9 else build_dir / f'reference-{spacing:g}.docx'
    make_reference_docx.build(target, spacing)
    return target


def build_docx(paper_dir, build_dir, refs, ID, spacing):
    raw = build_dir / f'{ID}_pandoc.docx'
    cmd = pandoc_common(paper_dir, build_dir, refs) + [
        '--to', 'docx', '--reference-doc', str(reference_doc(build_dir, spacing)), '-o', str(raw)]
    r = run(cmd, cwd=paper_dir, check=False)
    (build_dir / 'docx_pandoc.log').write_text(r.stderr)
    if r.returncode != 0:
        raise BuildError(f'pandoc (docx) failed:\n{r.stderr}')
    check_messages(r.stderr, 'docx')
    out = paper_dir / f'{out_stem(ID)}.docx'
    postprocess_docx(raw, out)
    return out


# ---- DOCX post-processing --------------------------------------------------------------------------
W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'


def _w(tag):
    return f'{{{W}}}{tag}'


SECTPR_ORDER = ['headerReference', 'footerReference', 'footnotePr', 'endnotePr', 'type', 'pgSz', 'pgMar',
                'paperSrc', 'pgBorders', 'lnNumType', 'pgNumType', 'cols', 'formProt', 'vAlign', 'noEndnote',
                'titlePg', 'textDirection', 'bidi', 'rtlGutter', 'docGrid', 'printerSettings', 'sectPrChange']
TBLPR_ORDER = ['tblStyle', 'tblpPr', 'tblOverlap', 'bidiVisual', 'tblStyleRowBandSize', 'tblStyleColBandSize',
               'tblW', 'jc', 'tblCellSpacing', 'tblInd', 'tblBorders', 'shd', 'tblLayout', 'tblCellMar', 'tblLook',
               'tblCaption', 'tblDescription']
TABLE_CELL_MARGIN = 40                 # twips, left and right, in data tables of the Word file
TABLE_TEXT_HALF_PT = b'18'             # size of the style "Table Text" in half-points (9 pt)
TCPR_ORDER = ['cnfStyle', 'tcW', 'gridSpan', 'hMerge', 'vMerge', 'tcBorders', 'shd', 'noWrap', 'tcMar',
              'textDirection', 'tcFitText', 'vAlign', 'hideMark']


def _put(parent, child, order):
    """Insert child into parent respecting the schema order (replacing an element with the same tag)."""
    local = child.tag.split('}')[1]
    for old in parent.findall(child.tag):
        parent.remove(old)
    idx = order.index(local)
    pos = len(parent)
    for i, el in enumerate(parent):
        name = el.tag.split('}')[1] if '}' in el.tag else el.tag
        if name in order and order.index(name) > idx:
            pos = i
            break
    parent.insert(pos, child)


def _el(tag, **attrs):
    from lxml import etree
    e = etree.Element(_w(tag))
    for k, v in attrs.items():
        e.set(_w(k), str(v))
    return e


def _borders(tag, spec):
    """tblBorders: every side explicit (nil where no rule). tcBorders: only the sides given, because a cell
    border (even nil) overrides the table border on that edge."""
    b = _el(tag)
    for side in ('top', 'left', 'bottom', 'right', 'insideH', 'insideV'):
        if tag == 'tcBorders' and side not in spec:
            continue
        v = spec.get(side)
        if v:
            b.append(_el(side, val='single', sz=v, space=0, color='auto'))
        else:
            b.append(_el(side, val='nil'))
    return b


CUSTOM_PROPS_EMPTY = (b'<?xml version="1.0" encoding="UTF-8"?><Properties '
                      b'xmlns="http://schemas.openxmlformats.org/officeDocument/2006/custom-properties" '
                      b'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes"/>')


def postprocess_docx(src, dst):
    from lxml import etree
    with zipfile.ZipFile(src) as z:
        items = [(i, z.read(i.filename)) for i in z.infolist()]
    doc_xml = next(data for i, data in items if i.filename == 'word/document.xml')
    root = etree.fromstring(doc_xml)
    body = root.find(_w('body'))

    # continuous line numbering in every section
    for sp in root.iter(_w('sectPr')):
        _put(sp, _el('lnNumType', countBy=1, distance=284, restart='continuous'), SECTPR_ORDER)

    # pandoc lays tables out for a 7920-twip (5.5 in) text block; rescale the grid to the real text width
    sect = body.find(_w('sectPr'))
    text_w = None
    if sect is not None and sect.find(_w('pgSz')) is not None and sect.find(_w('pgMar')) is not None:
        pg, mar = sect.find(_w('pgSz')), sect.find(_w('pgMar'))
        text_w = int(pg.get(_w('w'))) - int(mar.get(_w('left'))) - int(mar.get(_w('right')))

    n_eq = n_tab = 0
    for tbl in root.iter(_w('tbl')):
        grid = tbl.find(_w('tblGrid'))
        tblw = tbl.find(_w('tblPr') + '/' + _w('tblW')) if tbl.find(_w('tblPr')) is not None else None
        if text_w and grid is not None and tblw is not None and tblw.get(_w('type')) == 'pct':
            cols = grid.findall(_w('gridCol'))
            total = sum(int(c.get(_w('w'))) for c in cols) or 1
            frac = int(tblw.get(_w('w'))) / 5000.0
            for c in cols:
                c.set(_w('w'), str(int(round(int(c.get(_w('w'))) * text_w * frac / total))))
        tblpr = tbl.find(_w('tblPr'))
        if tblpr is None:
            tblpr = _el('tblPr')
            tbl.insert(0, tblpr)
        styles = {p.get(_w('val')) for p in tbl.iter(_w('pStyle'))}
        if 'EquationNumber' in styles:
            # numbered display equation: invisible 3-column table, equation centred, number flush right
            n_eq += 1
            _put(tblpr, _el('tblW', w=5000, type='pct'), TBLPR_ORDER)
            _put(tblpr, _borders('tblBorders', {}), TBLPR_ORDER)
            _put(tblpr, _el('tblLayout', type='fixed'), TBLPR_ORDER)
            _put(tblpr, _el('tblLook', firstRow=0, lastRow=0, firstColumn=0, lastColumn=0, noHBand=1, noVBand=1,
                            val='0000'), TBLPR_ORDER)
            for tc in tbl.iter(_w('tc')):
                tcpr = tc.find(_w('tcPr'))
                if tcpr is None:
                    tcpr = _el('tcPr')
                    tc.insert(0, tcpr)
                _put(tcpr, _el('vAlign', val='center'), TCPR_ORDER)
            continue
        n_tab += 1
        # data table: rules above and below the table and under the header row, nothing else
        _put(tblpr, _borders('tblBorders', {'top': 8, 'bottom': 8}), TBLPR_ORDER)
        _put(tblpr, _el('jc', val='center'), TBLPR_ORDER)
        # narrow cell margins (40 twips instead of Word's 108) and 9 pt table text (styles.xml, below): in the fixed
        # column grid a number must fit its cell, or Word breaks it over two lines (review of 5 Oct 2026)
        mar = _el('tblCellMar')
        for side in ('left', 'right'):
            mar.append(_el(side, w=TABLE_CELL_MARGIN, type='dxa'))
        _put(tblpr, mar, TBLPR_ORDER)
        rows = tbl.findall(_w('tr'))
        header_rows = [tr for tr in rows if tr.find(_w('trPr')) is not None
                       and tr.find(_w('trPr')).find(_w('tblHeader')) is not None]
        for tr in header_rows[-1:]:
            for tc in tr.findall(_w('tc')):
                tcpr = tc.find(_w('tcPr'))
                if tcpr is None:
                    tcpr = _el('tcPr')
                    tc.insert(0, tcpr)
                _put(tcpr, _borders('tcBorders', {'bottom': 4}), TCPR_ORDER)
                _put(tcpr, _el('vAlign', val='bottom'), TCPR_ORDER)
        for tr in rows:
            trpr = tr.find(_w('trPr'))
            if trpr is None:
                trpr = _el('trPr')
                tr.insert(0 if tr.find(_w('tblPrEx')) is None else 1, trpr)
            if trpr.find(_w('cantSplit')) is None:
                trpr.insert(0, _el('cantSplit'))
        for p in tbl.iter(_w('p')):
            ppr = p.find(_w('pPr'))
            if ppr is None:
                ppr = _el('pPr')
                p.insert(0, ppr)
            ps = ppr.find(_w('pStyle'))
            if ps is None:
                ps = _el('pStyle')
                ppr.insert(0, ps)
            ps.set(_w('val'), 'TableText')

    new_doc = etree.tostring(root, xml_declaration=True, encoding='UTF-8', standalone=True)
    tmp = dst.with_suffix('.tmp')
    n_props = 0
    with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as z:
        for info, data in items:
            if info.filename == 'word/document.xml':
                data = new_doc
            elif info.filename == 'word/styles.xml':
                data, n_tt = re.subn(rb'(<w:style [^>]*w:styleId="TableText".*?<w:sz w:val=")\d+(" ?/><w:szCs w:val=")\d+(")',
                                     lambda m: m.group(1) + TABLE_TEXT_HALF_PT + m.group(2) + TABLE_TEXT_HALF_PT + m.group(3),
                                     data, count=1, flags=re.S)
                if n_tt != 1:
                    raise BuildError('docx post-processing: the style "Table Text" was not found in styles.xml')
            elif info.filename == 'docProps/custom.xml':
                # pandoc stores the build metadata (crossref-map and csl paths of this machine, link-bibliography,
                # pdf-max-figure-height, xref-style) as custom document properties; the submitted Word file carries
                # none (ACS integration review, 30 Sep 2026)
                n_props = data.count(b'<property ')
                data = CUSTOM_PROPS_EMPTY
            z.writestr(info, data)
    tmp.replace(dst)
    say(f'  docx post-processing: {n_tab} table(s) ruled, {n_eq} numbered equation(s), line numbers on, '
        f'{n_props} build propert{"y" if n_props == 1 else "ies"} removed from docProps/custom.xml')


# ---- PDF --------------------------------------------------------------------------------------------
def _kpse(name, env):
    r = subprocess.run(['kpsewhich', name], capture_output=True, text=True, env=env)
    return r.stdout.strip() or None


def tex_env():
    env = dict(os.environ)
    if not _kpse('luaotfload-main.lua', env):
        local = REPO / 'tools' / 'texlive-extracted' / 'usr' / 'share' / 'texlive' / 'texmf-dist'
        if not local.exists():
            raise BuildError('luaotfload is not installed and tools/texlive-extracted is missing: fetch it with\n'
                             '  apt-get download texlive-luatex && dpkg-deb -x texlive-luatex_*.deb '
                             'tools/texlive-extracted')
        env['TEXMFHOME'] = str(local)
        if not _kpse('luaotfload-main.lua', env):
            raise BuildError('luaotfload still not found with TEXMFHOME=' + str(local))
    return env


def font_setup(env):
    """Paths for TeX Gyre Termes, the math font and the fallback chains (fonts are not installed system-wide)."""
    local_gyre = REPO / 'tools' / 'fonts' / 'extracted' / 'usr' / 'share' / 'texmf' / 'fonts' / 'opentype' / \
        'public' / 'tex-gyre'
    termes = local_gyre / 'texgyretermes-regular.otf'
    if not termes.exists():
        k = _kpse('texgyretermes-regular.otf', env)
        if not k:
            raise BuildError('TeX Gyre Termes not found (tools/fonts/extracted or TeX Live)')
        termes = Path(k)
    math = None
    for cand in (local_gyre / 'texgyretermes-math.otf', REPO / 'tools' / 'fonts' / 'texgyretermes-math.otf'):
        if cand.exists():
            math = cand
            break
    if math is None:
        k = _kpse('texgyretermes-math.otf', env) or _kpse('latinmodern-math.otf', env)
        if not k:
            raise BuildError('no OpenType math font (TeX Gyre Termes Math or Latin Modern Math)')
        math = Path(k)
    lib = REPO / 'tools' / 'fonts' / 'extracted' / 'usr' / 'share' / 'fonts' / 'truetype' / 'liberation'
    dejavu = Path('/usr/share/fonts/truetype/dejavu')
    chains = {}
    for shape, lib_name, dv_name in (('regular', 'LiberationSerif-Regular.ttf', 'DejaVuSerif.ttf'),
                                     ('italic', 'LiberationSerif-Italic.ttf', 'DejaVuSerif-Italic.ttf'),
                                     ('bold', 'LiberationSerif-Bold.ttf', 'DejaVuSerif-Bold.ttf'),
                                     ('bolditalic', 'LiberationSerif-BoldItalic.ttf', 'DejaVuSerif-BoldItalic.ttf')):
        c = []
        for f in (lib / lib_name, dejavu / dv_name, dejavu / 'DejaVuSerif.ttf', math):
            if f.exists() and f'[{f}]:mode=node;' not in c:
                c.append(f'[{f}]:mode=node;')
        chains[shape] = c
    return termes.parent, math, chains


def build_pdf(paper_dir, build_dir, refs, ID, spacing):
    env = tex_env()
    termes_dir, math, chains = font_setup(env)
    sp = {1.0: 'single', 1.5: 'onehalf', 2.0: 'double'}.get(round(spacing, 2))
    if sp is None:
        raise BuildError('--spacing must be 1, 1.5 or 2')
    var = ['--template', str(TEMPLATE), '--pdf-engine', 'lualatex',
           '-V', f'termes-dir={termes_dir}/', '-V', f'math-font-dir={math.parent}/', '-V', f'math-font-file={math.name}',
           '-V', f'spacing={sp}',
           '-V', f'figure-name={"Figure" if XREF == "acs" else "Fig."}',
           '--metadata', 'link-citations=true']
    if ID.endswith('_SI') and SI_STEM == '_Supplementary_Material':
        # Elsevier publishes supplementary files exactly as received: no line numbers in that PDF (the Word file of
        # the supplementary material, which is not uploaded, keeps the numbering of the template)
        var += ['-V', 'no-linenumbers=true']
    for shape, c in chains.items():
        for spec in c:
            var += ['-V', f'fallback-{shape}={spec}']
    common = pandoc_common(paper_dir, build_dir, refs)
    tex = build_dir / f'{out_stem(ID)}.tex'
    r = run(common + ['--to', 'latex', '--standalone'] + var + ['-o', str(tex)], cwd=paper_dir, env=env,
            check=False)
    if r.returncode != 0:
        raise BuildError(f'pandoc (latex) failed:\n{r.stderr}')
    check_messages(r.stderr, 'latex')
    out = paper_dir / f'{out_stem(ID)}.pdf'
    r = run(common + ['--to', 'pdf', '--verbose'] + var + ['-o', str(out)], cwd=paper_dir, env=env, check=False)
    log = build_dir / 'pdf_pandoc.log'
    log.write_text(r.stderr)
    missing = sorted(set(re.findall(r'Missing character: There is no (.+?) in font', r.stderr)))
    if missing:
        raise BuildError('Missing character(s) in the PDF: ' + ', '.join(missing) + f' (see {log})')
    if r.returncode != 0:
        tail = '\n'.join(l for l in r.stderr.splitlines() if l.startswith('!') or 'Error' in l)[-3000:]
        raise BuildError(f'pandoc/lualatex failed (full log: {log}):\n{tail}')
    # pandoc --verbose logs every LaTeX run; warnings that matter are those left after the final run
    runs = r.stderr.split('[makePDF] LaTeX run number')
    last = runs[-1] if len(runs) > 1 else r.stderr
    say(f'  [pdf] {max(1, len(runs) - 1)} LaTeX run(s); log in {rel(log)}')
    for l in r.stderr.splitlines():
        if l.startswith('[WARNING]'):
            say('  [pdf] ' + l)
    undefined = re.findall(r'LaTeX Warning: (?:Reference|Citation) [^\n]*undefined', last)
    if undefined:
        raise BuildError('undefined LaTeX reference(s) after the final run: ' + '; '.join(undefined[:5]))
    overfull = [float(x) for x in re.findall(r'Overfull \\hbox \(([\d.]+)pt too wide\)', last)]
    big = [x for x in overfull if x > 2.0]
    if big:
        say(f'  [pdf] {len(big)} line(s) stick into the margin by more than 2 pt (max {max(big):.1f} pt)')
    return out


# -------------------------------------------------------------------------------------------------
def main():
    global CSL, XREF, SI_STEM, ABBREV_STYLE
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paper')
    ap.add_argument('--docx', action='store_true', help='build the Word file')
    ap.add_argument('--pdf', action='store_true', help='build the PDF')
    ap.add_argument('--spacing', type=float, default=1.5, help='line spacing: 1.5 (default) or 2 (submission)')
    ap.add_argument('--no-bib', action='store_true', help='skip lit.py bib (use the existing references.json)')
    ap.add_argument('--journal', choices=sorted(JOURNALS), default=DEFAULT_JOURNAL,
                    help=f'citation profile (default {DEFAULT_JOURNAL})')
    a = ap.parse_args()
    CSL = JOURNALS[a.journal]['csl']
    XREF = JOURNALS[a.journal]['xref']
    SI_STEM = JOURNALS[a.journal]['si']
    ABBREV_STYLE = JOURNALS[a.journal]['abbrev']
    cite_mode = JOURNALS[a.journal]['citations']
    if not (a.docx or a.pdf):
        a.docx = a.pdf = True
    try:
        import lxml  # noqa: F401
    except ImportError:
        if PAPER_ENV_PY.exists() and Path(sys.executable).resolve() != PAPER_ENV_PY.resolve():
            os.execv(str(PAPER_ENV_PY), [str(PAPER_ENV_PY), __file__] + sys.argv[1:])
        say('lxml is needed (run with tools/paper-env/bin/python)')
        return 2
    try:
        paper_dir, ID = resolve(a.paper)
        if not (paper_dir / 'manuscript.md').exists():
            raise BuildError(f'{paper_dir / "manuscript.md"} does not exist')
        build_dir = paper_dir / 'build'
        build_dir.mkdir(exist_ok=True)
        say(f'== {ID}: {rel(paper_dir)}/manuscript.md  (journal profile {a.journal}: {JOURNALS[a.journal]["name"]}, '
            f'{rel(CSL)})')
        if not CSL.exists():
            raise BuildError(f'citation style {CSL} missing')
        if not a.no_bib:
            bibliography(paper_dir)
        refs = normalise_refs(paper_dir, build_dir)
        outputs = []
        if a.docx:
            say('-- DOCX')
            d = build_docx(paper_dir, build_dir, refs, ID, a.spacing)
            outputs.append(d)
            sys.path.insert(0, str(HERE))
            import docx_check
            rep = docx_check.check(d, crossref_map=build_dir / 'crossref.json', citations=cite_mode)
            say(docx_check.format_report(rep, indent='  '))
            if rep['errors']:
                raise BuildError(f'docx_check found {len(rep["errors"])} problem(s)')
        if a.pdf:
            say('-- PDF')
            outputs.append(build_pdf(paper_dir, build_dir, refs, ID, a.spacing))
        for o in outputs:
            say(f'OK  {rel(o)}  ({o.stat().st_size / 1024:.0f} kB)')
        return 0
    except BuildError as e:
        say(f'BUILD FAILED: {e}')
        return 1


if __name__ == '__main__':
    sys.exit(main())
