#!/usr/bin/env python3
"""Check a built manuscript .docx (standard library only).

  python3 paper/_shared/tools/docx_check.py <file.docx> [--crossref build/crossref.json] [--roundtrip] [--strict]
                                             [--citations auto|superscript|brackets]

Reports
  * figures (drawings, "Fig. N." captions), tables (excluding the invisible equation tables), numbered and
    unnumbered display equations, inline math
  * cross-references in the text (Fig./Table/Eq. numbers) against what exists; figures/tables never referenced
  * leftovers of unresolved markup: '??', '@fig:', '@tbl:', '@eq:', '@sec:', '[@', '{#'
  * citations: numbers used in the text vs the reference list (every reference cited, every citation listed,
    first appearances in order 1, 2, 3, ...). Two forms (--citations, default auto = from the document):
      superscript  ACS profile: citation numbers are the runs in the character style "Citation Number" that
                   build.py/bibfix.lua puts on every citation ("salt.^1,3-5^"); they must be superscript and read
                   n, n-m, n,m,...; reference list entries start "(n)". Other superscripts (cm⁻¹, R²₂(8)) are
                   ignored, and citation runs are left out of the text used for the Fig./Table/Eq. checks, so
                   that "Table S17^5,6^" is not read as Table S175.
      brackets     Elsevier profile: "[1]", "[2,3]", "[4-6]" in the text; entries start "[n]".
  * images: pixel size, placed size and effective resolution (warning below 300 dpi)
  * line numbering and a page-number field in the footer
  * --roundtrip: converts the .docx back to Markdown with pandoc and checks that captions, citations and the
    reference list survive (superscript form: read with docx+styles; the styled citations that come back must be
    the same as in the document, number for number, and cover 1..N. Their order is reported but not checked, because
    pandoc's Markdown puts a table caption after its table; the order is checked on the document itself.)
Exit status 1 when an error is found (with --strict also on warnings).
"""
import argparse
import json
import re
import struct
import subprocess
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

NS = {
    'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'm': 'http://schemas.openxmlformats.org/officeDocument/2006/math',
    'wp': 'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
    'a': 'http://schemas.openxmlformats.org/drawingml/2006/main',
    'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
    'pr': 'http://schemas.openxmlformats.org/package/2006/relationships',
}
W = '{%s}' % NS['w']
M = '{%s}' % NS['m']
EMU_PER_CM = 360000
REPO = Path(__file__).resolve().parents[3]

CITE_RE = re.compile(r'\[(\d+(?:\s*[\u2013\-]\s*\d+)?(?:\s*,\s*\d+(?:\s*[\u2013\-]\s*\d+)?)*)\]')
SUP_CITE_RE = re.compile(r'^(\d+(?:[\u2013\-]\d+)?(?:,\d+(?:[\u2013\-]\d+)?)*)$')
CITE_STYLE_IDS = {'CitationNumber'}          # pandoc's style id for custom-style="Citation Number"
CITE_STYLE_NAME = 'Citation Number'
NO_XREF_STYLES = {'ContentsList'}            # SI contents list (crossref.lua): not a cross-reference in the text
XREF_RE = {
    'fig': re.compile(r'\b(?:Figs?\.|Figures?)\u00a0?\s?((?:S?\d+[a-z]?(?:\u2013S?\d+[a-z]?)?(?:,\s|\sand\s)?)+)'),
    'tbl': re.compile(r'\bTables?\u00a0?\s?((?:S?\d+(?:\u2013S?\d+)?(?:,\s|\sand\s)?)+)'),
    # "Eq. (1)" (Elsevier profile) or "eq 1" (ACS profile, crossref.lua xref-style=acs)
    'eq': re.compile(r'\b(?:Eqs?\.\u00a0?\s?((?:\(S?\d+\)(?:\u2013\(S?\d+\))?(?:,\s|\sand\s)?)+)|'
                     r'eqs?\u00a0((?:S?\d+(?:\u2013S?\d+)?(?:,\s|\sand\s)?)+))'),
}
LEFTOVER_RE = re.compile(r'\?\?|@(?:fig|tbl|eq|sec):|\[@|\{#')


def expand(spec):
    out = []
    for part in re.split(r'\s*,\s*', spec.strip()):
        m = re.match(r'(\d+)\s*[\u2013\-]\s*(\d+)$', part)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            out.extend(range(a, b + 1) if b >= a else [a, b])
        elif part.isdigit():
            out.append(int(part))
    return out


def xref_numbers(s, prefix=None):
    """Numbers referred to in a 'Fig./Table' match. prefix '' (main text) or 'S' (Supporting Information):
    only references in the document's own series are returned, so that a main text may cite 'Fig. S3' and an
    SI may cite the main text's 'Fig. 3' without either being counted as a figure of its own. None = all."""
    nums = []
    for tok in re.split(r',\s|\sand\s', s):
        tok = tok.strip().strip('()')
        if prefix is not None and tok.lstrip('(').startswith('S') != (prefix == 'S'):
            continue
        m = re.match(r'\(?S?(\d+)[a-z]?\)?(?:\u2013\(?S?(\d+)[a-z]?\)?)?$', tok)
        if m:
            a = int(m.group(1))
            b = int(m.group(2)) if m.group(2) else a
            nums.extend(range(a, b + 1))
    return nums


def image_size(data):
    """(width, height) in pixels for PNG / JPEG / GIF / TIFF-less data, or None."""
    if data[:8] == b'\x89PNG\r\n\x1a\n':
        return struct.unpack('>II', data[16:24])
    if data[:6] in (b'GIF87a', b'GIF89a'):
        return struct.unpack('<HH', data[6:10])
    if data[:2] == b'\xff\xd8':
        i = 2
        while i < len(data) - 9:
            if data[i] != 0xFF:
                i += 1
                continue
            marker = data[i + 1]
            if marker in (0xC0, 0xC1, 0xC2, 0xC3, 0xC5, 0xC6, 0xC7, 0xC9, 0xCA, 0xCB, 0xCD, 0xCE, 0xCF):
                h, w = struct.unpack('>HH', data[i + 5:i + 9])
                return w, h
            seg = struct.unpack('>H', data[i + 2:i + 4])[0]
            i += 2 + seg
    return None


def run_style(r):
    rs = r.find('w:rPr/w:rStyle', NS)
    return rs.get(W + 'val') if rs is not None else ''


def is_cite_run(r):
    return run_style(r) in CITE_STYLE_IDS


def run_is_superscript(r, style_sup):
    va = r.find('w:rPr/w:vertAlign', NS)
    if va is not None:
        return va.get(W + 'val') == 'superscript'
    return style_sup


def para_text(p, skip_cites=False):
    """Visible text of a paragraph (runs, including hyperlinks), math rendered as its m:t text. skip_cites: leave
    out the runs in the "Citation Number" style (superscript citations)."""
    parts = []
    skip = set()
    if skip_cites:
        for r in p.iter(W + 'r'):
            if is_cite_run(r):
                skip.update(id(x) for x in r.iter())
    for el in p.iter():
        if id(el) in skip:
            continue
        if el.tag == W + 't' or el.tag == M + 't':
            parts.append(el.text or '')
        elif el.tag == W + 'tab':
            parts.append('\t')
        elif el.tag in (W + 'br', W + 'cr'):
            parts.append('\n')
        elif el.tag == W + 'noBreakHyphen':
            parts.append('-')
    return ''.join(parts)


def pstyle(p):
    ps = p.find('w:pPr/w:pStyle', NS)
    return ps.get(W + 'val') if ps is not None else ''


def sup_citations(p, style_sup, errors):
    """Citation strings of a paragraph in the superscript form: consecutive runs in the citation style form one
    citation. A citation run that is not superscript is an error."""
    groups, cur = [], None
    for r in p.iter(W + 'r'):
        if is_cite_run(r):
            txt = ''.join(t.text or '' for t in r.iter(W + 't'))
            if not run_is_superscript(r, style_sup):
                errors.append(f'citation {txt!r} is not superscript')
            if cur is None:
                cur = []
                groups.append(cur)
            cur.append(txt)
        elif ''.join(t.text or '' for t in r.iter(W + 't')) != '' or r.find('w:drawing', NS) is not None:
            cur = None
    return [''.join(g) for g in groups]


def citation_mode(doc, styles_xml, bib):
    """'superscript' when the document has runs in the citation style or '(n)' reference entries, else 'brackets'"""
    if any(is_cite_run(r) for r in doc.iter(W + 'r')):
        return 'superscript'
    if bib and re.match(r'\s*\(\d+\)', bib[0]):
        return 'superscript'
    return 'brackets'


def style_is_superscript(styles_xml):
    """True when the "Citation Number" character style itself is superscript (runs may inherit it)."""
    if not styles_xml:
        return False
    st = ET.fromstring(styles_xml)
    for s in st.findall('w:style', NS):
        if s.get(W + 'styleId') in CITE_STYLE_IDS:
            va = s.find('w:rPr/w:vertAlign', NS)
            return va is not None and va.get(W + 'val') == 'superscript'
    return False


def check(path, crossref_map=None, roundtrip=False, citations='auto'):
    path = Path(path)
    errors, warnings = [], []
    rep = {'file': str(path), 'errors': errors, 'warnings': warnings}
    with zipfile.ZipFile(path) as z:
        names = set(z.namelist())
        doc = ET.fromstring(z.read('word/document.xml'))
        styles_xml = z.read('word/styles.xml') if 'word/styles.xml' in names else b''
        rels = ET.fromstring(z.read('word/_rels/document.xml.rels'))
        relmap = {r.get('Id'): r.get('Target') for r in rels.findall('pr:Relationship', NS)}
        media = {n: z.read(n) for n in names if n.startswith('word/media/')}
        footers = [z.read(n).decode('utf-8', 'ignore') for n in names if re.match(r'word/footer\d*\.xml$', n)]
    body = doc.find('w:body', NS)

    # ---- tables: data tables vs equation tables
    eq_tables, data_tables = [], []
    for t in doc.iter(W + 'tbl'):
        styles = {pstyle(p) for p in t.iter(W + 'p')}
        (eq_tables if 'EquationNumber' in styles else data_tables).append(t)
    in_eq_table = {id(p) for t in eq_tables for p in t.iter(W + 'p')}
    in_table = {id(p) for t in data_tables for p in t.iter(W + 'p')}

    paras = list(doc.iter(W + 'p'))
    bib_all = [para_text(p) for p in paras if pstyle(p) == 'Bibliography']
    mode = citation_mode(doc, styles_xml, bib_all) if citations == 'auto' else citations
    if mode not in ('superscript', 'brackets'):
        raise ValueError(f'unknown citation form {citations!r}')
    sup = mode == 'superscript'
    style_sup = style_is_superscript(styles_xml)
    bib, body_paras, fig_caps, tab_caps = [], [], [], []
    para_cites = []                                  # citation strings per body paragraph, in document order
    for p in paras:
        st = pstyle(p)
        txt = para_text(p, skip_cites=sup)
        if st == 'Bibliography':
            bib.append(txt)
        elif st == 'ImageCaption':
            fig_caps.append(txt)
        elif st == 'TableCaption':
            tab_caps.append(txt)
        if st != 'Bibliography' and id(p) not in in_eq_table:
            body_paras.append((st, txt, id(p) in in_table))
            if sup:
                para_cites.append(sup_citations(p, style_sup, errors))
            else:
                para_cites.append([m.group(1) for m in CITE_RE.finditer(txt)])

    # ---- figures
    drawings = list(doc.iter('{%s}inline' % NS['wp'])) + list(doc.iter('{%s}anchor' % NS['wp']))
    fig_nums = []
    for c in fig_caps:
        m = re.match(r'(?:Fig\.|Figure)\u00a0?\s?S?(\d+)\.', c)
        if m:
            fig_nums.append(int(m.group(1)))
        else:
            errors.append(f'figure caption without "Fig. N." or "Figure N.": {c[:60]!r}')
    if fig_nums != list(range(1, len(fig_nums) + 1)):
        errors.append(f'figure captions not numbered 1..N in order: {fig_nums}')
    tab_nums = []
    for c in tab_caps:
        m = re.match(r'Table\u00a0?\s?S?(\d+)\.', c)
        if m:
            tab_nums.append(int(m.group(1)))
        else:
            errors.append(f'table caption without "Table N.": {c[:60]!r}')
    if tab_nums != list(range(1, len(tab_nums) + 1)):
        errors.append(f'table captions not numbered 1..N in order: {tab_nums}')
    if len(tab_caps) != len(data_tables):
        warnings.append(f'{len(data_tables)} tables but {len(tab_caps)} table captions')

    # ---- equations
    eq_nums = []
    for t in eq_tables:
        for p in t.iter(W + 'p'):
            if pstyle(p) == 'EquationNumber':
                m = re.match(r'\(S?(\d+)\)$', para_text(p).strip())
                if m:
                    eq_nums.append(int(m.group(1)))
    if eq_nums != list(range(1, len(eq_nums) + 1)):
        errors.append(f'equation numbers not 1..N in order: {eq_nums}')
    omath_para = list(doc.iter(M + 'oMathPara'))
    inline_math = [m for m in doc.iter(M + 'oMath')]
    n_display = len(omath_para)
    n_inline = len(inline_math) - sum(len(list(mp.iter(M + 'oMath'))) for mp in omath_para)

    rep['figures'] = {'drawings': len(drawings), 'captions': fig_nums}
    rep['tables'] = {'tables': len(data_tables), 'captions': tab_nums}
    rep['equations'] = {'numbered': eq_nums, 'display_math': n_display, 'inline_math': n_inline}

    # ---- cross-references in the running text (captions excluded)
    running = ' '.join(txt for st, txt, intab in body_paras
                       if st not in ('ImageCaption', 'TableCaption') and st not in NO_XREF_STYLES and not intab)
    used = {}
    own_prefix = 'S' if any(re.match(r'(?:Fig\.|Figure|Table)\u00a0?\s?S\d', c) for c in fig_caps + tab_caps) else ''
    for kind, rx in XREF_RE.items():
        nums = set()
        for m in rx.finditer(running):
            spec = next(g for g in m.groups() if g is not None)
            nums.update(xref_numbers(spec, own_prefix if kind != 'eq' else None))
        used[kind] = sorted(nums)
    rep['xrefs'] = used
    for kind, have in (('fig', len(fig_nums)), ('tbl', len(tab_nums)), ('eq', len(eq_nums))):
        name = {'fig': 'Fig.', 'tbl': 'Table', 'eq': 'Eq.'}[kind]
        bad = [n for n in used[kind] if n > have or n < 1]
        if bad:
            errors.append(f'text refers to {name} {bad} but only {have} exist')
        unref = [n for n in range(1, have + 1) if n not in used[kind]]
        if unref and kind != 'eq':
            warnings.append(f'{name} {unref} never referred to in the text')

    # ---- leftovers
    for st, txt, _ in body_paras:
        for m in LEFTOVER_RE.finditer(txt):
            errors.append(f'unresolved markup {m.group(0)!r} in: {txt[max(0, m.start() - 40):m.end() + 40]!r}')
    for txt in bib:
        if '??' in txt:
            errors.append(f'"??" in reference list entry {txt[:60]!r}')

    # ---- citations
    cited_order, cited = [], set()
    n_cites = 0
    all_cites = []
    for cites in para_cites:
        for c in cites:
            n_cites += 1
            all_cites.append(c)
            if sup and not SUP_CITE_RE.match(c):
                errors.append(f'malformed superscript citation {c!r} (expected n, n\u2013m or n,m,...)')
                continue
            for n in expand(c):
                if n not in cited:
                    cited_order.append(n)
                cited.add(n)
    if sup:
        stray = [m.group(0) for st, txt, _ in body_paras for m in CITE_RE.finditer(txt)]
        if stray:
            warnings.append(f'{len(stray)} bracketed number(s) such as {stray[:3]} in a superscript-citation '
                            f'document (a citation that escaped the style?)')
    entry_re = r'\s*\((\d+)\)' if sup else r'\s*\[(\d+)\]'
    entries = []
    for txt in bib:
        m = re.match(entry_re, txt)
        if m:
            entries.append(int(m.group(1)))
        else:
            errors.append(f'reference list entry without {"(n)" if sup else "[n]"}: {txt[:60]!r}')
    rep['citations'] = {'form': mode, 'citations': n_cites, 'cited': sorted(cited),
                        'first_order_ok': cited_order == sorted(cited_order), 'references': len(entries)}
    if entries != list(range(1, len(entries) + 1)):
        errors.append(f'reference list not numbered 1..N: {entries[:10]}...')
    not_listed = sorted(n for n in cited if n not in entries)
    not_cited = sorted(n for n in entries if n not in cited)
    if not_listed:
        errors.append(f'cited but not in the reference list: {not_listed}')
    if not_cited:
        errors.append(f'in the reference list but never cited: {not_cited}')
    if cited_order != list(range(1, len(cited_order) + 1)):
        errors.append(f'citations do not first appear in the order 1, 2, 3, ...: {cited_order[:15]}')
    no_doi = [txt[:50] for txt in bib if 'doi.org/' not in txt]
    if no_doi:
        warnings.append(f'{len(no_doi)} reference(s) without a DOI: {no_doi[:3]}')

    # ---- images
    imgs = []
    for d in drawings:
        ext = d.find('wp:extent', NS)
        blip = d.find('.//a:blip', NS)
        if ext is None or blip is None:
            continue
        target = relmap.get(blip.get('{%s}embed' % NS['r']), '')
        data = media.get('word/' + target, b'')
        px = image_size(data)
        cx, cy = int(ext.get('cx')), int(ext.get('cy'))
        w_cm, h_cm = cx / EMU_PER_CM, cy / EMU_PER_CM
        info = {'media': target, 'px': px, 'cm': (round(w_cm, 2), round(h_cm, 2))}
        if px:
            dpi = px[0] / (w_cm / 2.54)
            info['dpi'] = round(dpi)
            if dpi < 300:
                warnings.append(f'{target}: {dpi:.0f} dpi at the placed size ({w_cm:.1f} cm wide); journals ask '
                                f'for >= 300 (photos) or >= 600 (line art)')
            if abs(px[1] / px[0] - h_cm / w_cm) > 0.01:
                errors.append(f'{target}: aspect ratio distorted')
        if w_cm > 16.0 + 0.01:
            warnings.append(f'{target}: placed {w_cm:.1f} cm wide (> 16 cm)')
        imgs.append(info)
    rep['images'] = imgs

    # ---- line numbers, footer page field
    ln = [s.find('w:lnNumType', NS) for s in doc.iter(W + 'sectPr')]
    rep['line_numbers'] = bool(ln) and all(x is not None and x.get(W + 'restart') == 'continuous' for x in ln)
    if not rep['line_numbers']:
        errors.append('line numbering (w:lnNumType restart=continuous) missing')
    rep['page_numbers'] = any('PAGE' in f for f in footers) and \
        any(s.find('w:footerReference', NS) is not None for s in doc.iter(W + 'sectPr'))
    if not rep['page_numbers']:
        errors.append('no PAGE field in a footer referenced by the section')

    # ---- crossref map from the build (optional)
    if crossref_map and Path(crossref_map).exists():
        cm = json.loads(Path(crossref_map).read_text())
        kinds = {}
        for lab, v in cm.items():
            kinds.setdefault(v['kind'], []).append(v['num'])
        for kind, have in (('fig', fig_nums), ('tbl', tab_nums), ('eq', eq_nums)):
            labelled = sorted(int(n) for n in kinds.get(kind, []))
            missing = [n for n in labelled if n not in have]
            if missing:
                errors.append(f'{kind} labels numbered {missing} by crossref.lua are absent from the docx')
        rep['labels'] = len(cm)

    # ---- round trip through pandoc
    if roundtrip:
        rep['roundtrip'] = roundtrip_check(path, len(fig_nums), len(tab_nums), len(entries), errors,
                                           sup_cites=all_cites if sup else None)

    words = sum(len(t.split()) for st, t, intab in body_paras if not intab)
    rep['words'] = words
    return rep


def roundtrip_check(path, n_fig, n_tab, n_ref, errors, sup_cites=None):
    """sup_cites: the superscript citations found in the document, in order (None for the bracket form)."""
    pandoc = REPO / 'tools' / 'bin' / 'pandoc'
    exe = str(pandoc) if pandoc.exists() else 'pandoc'
    r = subprocess.run([exe, str(path), '-f', 'docx', '-t', 'markdown', '--wrap=none'],
                       capture_output=True, text=True)
    if r.returncode != 0:
        errors.append('pandoc could not read the docx back: ' + r.stderr[:300])
        return {}
    md = r.stdout
    if sup_cites is None:
        refs = len(re.findall(r'^\[\d+\]', md, flags=re.M)) + len(re.findall(r'^\\\[\d+\\\]', md, flags=re.M))
        cites = len(re.findall(r'\\?\[\d+(?:[\u2013,\-]\d+)*\\?\]', md))
    else:
        refs = len(re.findall(r'^\\?\(\d+\\?\)\s', md, flags=re.M))
        cites = None
    out = {
        'fig_captions': len(re.findall(r'\*\*(?:Fig\.|Figure)\u00a0S?\d+\.\*\*|\*\*(?:Fig\.|Figure) S?\d+\.\*\*', md)),
        'table_captions': len(re.findall(r'\*\*Table\u00a0S?\d+\.\*\*|\*\*Table S?\d+\.\*\*', md)),
        'reference_entries': refs,
        'citations': cites,
        'math': md.count('$'),
        'length': len(md),
    }
    if sup_cites is not None:
        # read again keeping the character styles: every citation must come back as a superscript in the
        # "Citation Number" style, as many as in the document, first appearances 1..N in order
        r2 = subprocess.run([exe, str(path), '-f', 'docx+styles', '-t', 'markdown-smart', '--wrap=none'],
                            capture_output=True, text=True)
        if r2.returncode != 0:
            errors.append('pandoc could not read the docx back with its styles: ' + r2.stderr[:300])
        else:
            # [^1,2^]{custom-style="Citation Number"}, or raw HTML where Markdown cannot express the span
            # (e.g. directly after a closing parenthesis): <span data-custom-style="Citation Number"><sup>22</sup></span>
            name = re.escape(CITE_STYLE_NAME)
            rx = re.compile(r'\[\^([^\^\]]*)\^\]\{custom-style="' + name + r'"\}|'
                            r'<span data-custom-style="' + name + r'"><sup>([^<]*)</sup></span>')
            found = [re.sub(r'\\(.)', r'\1', m.group(1) if m.group(1) is not None else m.group(2))
                     for m in rx.finditer(r2.stdout)]
            n_styled = r2.stdout.count('custom-style="' + CITE_STYLE_NAME + '"')
            if n_styled != len(found):
                errors.append(f'round trip: {n_styled - len(found)} citation(s) in the "{CITE_STYLE_NAME}" style '
                              f'come back in an unexpected form')
            out['citations'] = len(found)
            order, seen = [], set()
            for c in found:
                if not SUP_CITE_RE.match(c):
                    errors.append(f'round trip: citation {c!r} is not a superscript number list')
                    continue
                for n in expand(c):
                    if n not in seen:
                        order.append(n)
                        seen.add(n)
            out['citations_first_order'] = 'yes' if order == list(range(1, len(order) + 1)) else \
                'no (informational: Markdown moves table captions after their tables; the order is checked on the document)'
            if len(found) != len(sup_cites):
                errors.append(f'round trip: {len(found)} superscript citations survive, expected {len(sup_cites)}')
            elif sorted(found) != sorted(sup_cites):
                diff = sorted(set(sup_cites) ^ set(found))
                errors.append(f'round trip: the citations that survive differ from those of the document: {diff[:10]}')
            if sorted(seen) != list(range(1, n_ref + 1)):
                errors.append(f'round trip: citations cover {compress(seen)}, the list has 1-{n_ref}')
    if out['fig_captions'] != n_fig:
        errors.append(f'round trip: {out["fig_captions"]} figure captions survive, expected {n_fig}')
    if out['table_captions'] != n_tab:
        errors.append(f'round trip: {out["table_captions"]} table captions survive, expected {n_tab}')
    if out['reference_entries'] < n_ref:
        errors.append(f'round trip: {out["reference_entries"]} reference entries survive, expected {n_ref}')
    return out


def format_report(rep, indent=''):
    L = []
    f, t, e, c = rep['figures'], rep['tables'], rep['equations'], rep['citations']
    L.append(f'docx check: {rep["file"]}  (~{rep["words"]} words outside tables)')
    L.append(f'  figures: {f["drawings"]} image(s); captions Fig. {f["captions"]}')
    L.append(f'  tables: {t["tables"]}; captions Table {t["captions"]}')
    L.append(f'  equations: numbered {e["numbered"]}; display math objects {e["display_math"]}; '
             f'inline math {e["inline_math"]}')
    x = rep['xrefs']
    L.append(f'  cross-references in text: Fig. {x["fig"]}  Table {x["tbl"]}  Eq. {x["eq"]}')
    L.append(f'  citations ({c.get("form", "brackets")}): {c.get("citations", "?")} citation(s) citing '
             f'{compress(c["cited"])}; reference list {c["references"]} entries; '
             f'first appearances in order: {"yes" if c["first_order_ok"] else "NO"}')
    for im in rep['images']:
        px = f'{im["px"][0]}x{im["px"][1]} px' if im.get('px') else '? px'
        L.append(f'  image {im["media"]}: {px}, placed {im["cm"][0]} x {im["cm"][1]} cm -> {im.get("dpi", "?")} dpi')
    L.append(f'  line numbers: {"on" if rep["line_numbers"] else "OFF"}; page numbers: '
             f'{"on" if rep["page_numbers"] else "OFF"}')
    if 'roundtrip' in rep:
        L.append(f'  round trip (docx -> markdown): {rep["roundtrip"]}')
    for w in rep['warnings']:
        L.append(f'  WARNING: {w}')
    for er in rep['errors']:
        L.append(f'  ERROR: {er}')
    L.append(f'  RESULT: {"OK" if not rep["errors"] else str(len(rep["errors"])) + " error(s)"}'
             f'{", " + str(len(rep["warnings"])) + " warning(s)" if rep["warnings"] else ""}')
    return '\n'.join(indent + l for l in L)


def compress(nums):
    nums = sorted(nums)
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(f'{nums[i]}' if i == j else f'{nums[i]}-{nums[j]}')
        i = j + 1
    return '[' + ', '.join(out) + ']'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('docx')
    ap.add_argument('--crossref', help='crossref.json written by the build (label -> number)')
    ap.add_argument('--roundtrip', action='store_true', help='also convert back to Markdown with pandoc')
    ap.add_argument('--strict', action='store_true', help='fail on warnings too')
    ap.add_argument('--json', action='store_true', help='print the report as JSON')
    ap.add_argument('--citations', choices=('auto', 'superscript', 'brackets'), default='auto',
                    help='citation form: superscript (ACS), brackets (Elsevier) or auto (from the document)')
    a = ap.parse_args()
    rep = check(a.docx, a.crossref, a.roundtrip, a.citations)
    print(json.dumps(rep, indent=1, ensure_ascii=False) if a.json else format_report(rep))
    return 1 if rep['errors'] or (a.strict and rep['warnings']) else 0


if __name__ == '__main__':
    sys.exit(main())
