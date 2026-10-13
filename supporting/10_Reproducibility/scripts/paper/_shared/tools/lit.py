#!/usr/bin/env python3
"""Literature helper for the SAPH manuscript and its Supporting Information.

Every reference that reaches the paper passes through this script: metadata is taken from Crossref, never
typed by hand, so a citation cannot carry an invented title, volume or page. Adapted from cqc-dft
(papers/_shared/tools/lit.py): one literature pool for this project (lit/candidates.jsonl) instead of one per
paper, a stricter retraction check, a must-cite flag, a batch add, a full re-verification and a Markdown listing.

  lit.py search "query" [--from 2020] [--to 2026] [--n 25] [--abstracts]   OpenAlex keyword search
  lit.py check DOI [DOI ...]                        Crossref lookup + retraction check; exit 1 if any fails
  lit.py abstract DOI                               abstract text (OpenAlex), for checking relevance
  lit.py add DOI --topic T --claim "..." [--relevance 1-5] [--section S] [--must-cite]
                                                    verify, then append to lit/candidates.jsonl
  lit.py add-many FILE.json                         the same for a list of {doi, topic, claim, relevance,
                                                    section, must_cite}; one line of output per DOI
  lit.py drop DOI --reason "..."                    mark a candidate as rejected (kept for the record)
  lit.py claim KEY --old "..." --new "..."          replace one claim of a candidate (its exact text), e.g. when
                                                    a claim first written from the title is rewritten from the
                                                    abstract; the Crossref metadata are not touched
  lit.py list [--all]                               candidates with key, year, topic, claim
  lit.py stats                                      counts by year, topic, journal; must-cite status
  lit.py verify [--no-refresh]                      re-fetch every active candidate from Crossref and OpenAlex
                                                    (bypassing the cache) and redo the retraction check; exit 1
                                                    if any DOI no longer resolves or is now flagged
  lit.py table                                      write lit/REFS.md (grouped listing for the authors)
  lit.py bib [--manuscript FILE] [--out FILE]       CSL-JSON (pandoc) of the keys cited in FILE
                                                    (default paper/main/manuscript.md -> paper/main/references.json);
                                                    applies lit/csl_overrides.json corrections (each needs "_why")

Retraction check (a DOI with any of these is refused by add/add-many and reported by check/verify):
  * OpenAlex is_retracted;
  * Crossref 'update-to' or 'updated-by' entries of type retraction, partial_retraction, withdrawal or removal
    (Crossref merges the Retraction Watch database into 'updated-by');
  * any Crossref notice returned by works?filter=updates:DOI with one of those types;
  * a title starting "Retracted" / "Retraction" / "Withdrawn".
An expression of concern is refused too unless --allow-concern is given. Corrections, errata and addenda are
not refusals; they are stored in the candidate's "updates" field so that the author can cite the correction.

Citation keys are deterministic (surnameYEARfirstword as in cqc-dft, but with markup such as <i> or <scp>
stripped from the title first), so they can be written into the manuscript before the bibliography is generated.
Network calls are spaced by a cross-process lock, and every Crossref/OpenAlex record is cached in lit/cache/
(git-ignored).
"""
import argparse, datetime, fcntl, json, re, sys, time, unicodedata
from pathlib import Path
import urllib.parse, urllib.request, urllib.error

REPO = Path(__file__).resolve().parents[3]           # repo root (paper/_shared/tools/lit.py)
LIT = REPO / 'lit'
CAND = LIT / 'candidates.jsonl'
CACHE = LIT / 'cache'
LOCK = CACHE / '.api.lock'
MAIN_MS = REPO / 'paper' / 'main' / 'manuscript.md'
# No e-mail address in the User-Agent: the polite pool is not needed at this volume.
UA = 'cqc-sfa-paper/1.0 (literature verification script)'
REFUSE_TYPES = ('retraction', 'partial_retraction', 'withdrawal', 'removal')
CONCERN_TYPES = ('expression_of_concern',)


class RateLimited(Exception):
    pass


def _get(url, accept='application/json', min_interval=0.35, tries=6, no_wait_429=False):
    """GET with a cross-process lock and spacing, so parallel agents do not trip the API rate limits."""
    CACHE.mkdir(parents=True, exist_ok=True)
    for attempt in range(tries):
        with open(LOCK, 'a+') as lk:
            fcntl.flock(lk, fcntl.LOCK_EX)
            lk.seek(0)
            try:
                last = float(lk.read().strip() or 0)
            except ValueError:
                last = 0.0
            wait = last + min_interval - time.time()
            if wait > 0:
                time.sleep(wait)
            lk.seek(0); lk.truncate(); lk.write(str(time.time())); lk.flush()
            fcntl.flock(lk, fcntl.LOCK_UN)
        req = urllib.request.Request(url, headers={'User-Agent': UA, 'Accept': accept})
        try:
            with urllib.request.urlopen(req, timeout=40) as r:
                return r.read().decode('utf-8')
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if e.code == 429 and no_wait_429:
                raise RateLimited(url)
            if e.code in (429, 500, 502, 503, 504):
                delay = 2 ** attempt
                try:
                    body = json.loads(e.read().decode('utf-8', 'ignore') or '{}')
                    delay = max(delay, min(float(body.get('retryAfter', 0)), 65))
                except Exception:
                    pass
                ra = e.headers.get('Retry-After') if e.headers else None
                if ra and ra.isdigit():
                    delay = max(delay, min(int(ra), 65))
                time.sleep(delay)
                continue
            raise
        except (urllib.error.URLError, TimeoutError):
            time.sleep(2 ** attempt)
    raise RuntimeError(f'giving up on {url}')


def norm_doi(doi):
    doi = doi.strip()
    doi = re.sub(r'^(https?://)?(dx\.)?doi\.org/', '', doi, flags=re.I)
    doi = re.sub(r'^doi:\s*', '', doi, flags=re.I)
    return doi.rstrip('.').lower()


def _safe(doi):
    return re.sub(r'[^a-z0-9._-]+', '_', doi)


def _cached(kind, doi, url, refresh=False, accept='application/json', min_interval=0.35, unwrap=None):
    """Fetch url once and keep it in lit/cache/<kind>/<doi>.json; refresh=True re-fetches."""
    f = CACHE / kind / (_safe(norm_doi(doi)) + '.json')
    if f.exists() and not refresh:
        return json.loads(f.read_text())
    txt = _get(url, accept=accept, min_interval=min_interval)
    if txt is None:
        return None
    d = json.loads(txt)
    if unwrap:
        d = unwrap(d)
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(d))
    return d


def crossref(doi, refresh=False):
    """Crossref record (message) for a DOI, cached. None if Crossref does not know it."""
    return _cached('crossref', doi, 'https://api.crossref.org/works/' + urllib.parse.quote(norm_doi(doi), safe=''),
                   refresh, unwrap=lambda d: d['message'])


def csl(doi, refresh=False):
    """CSL-JSON straight from Crossref's transform endpoint, cached."""
    return _cached('csl', doi, 'https://api.crossref.org/works/' + urllib.parse.quote(norm_doi(doi), safe='') +
                   '/transform/application/vnd.citationstyles.csl+json', refresh,
                   accept='application/vnd.citationstyles.csl+json')


def openalex(doi, refresh=False):
    return _cached('openalex', doi, 'https://api.openalex.org/works/doi:' + urllib.parse.quote(norm_doi(doi), safe='/'),
                   refresh, min_interval=0.15)


def update_notices(doi, refresh=False):
    """Crossref works that declare they update this DOI (retraction notices, errata, ...), cached."""
    q = urllib.parse.urlencode({'filter': f'updates:{norm_doi(doi)}', 'select': 'DOI,type,update-to,title',
                                'rows': '50'})
    d = _cached('updates', doi, 'https://api.crossref.org/works?' + q, refresh,
                unwrap=lambda d: d['message'].get('items', []))
    return d or []


def abstract_text(oa):
    inv = (oa or {}).get('abstract_inverted_index')
    if not inv:
        return ''
    pos = {}
    for w, idx in inv.items():
        for i in idx:
            pos[i] = w
    return ' '.join(pos[i] for i in sorted(pos))


def _ascii(s):
    return unicodedata.normalize('NFKD', s).encode('ascii', 'ignore').decode()


STOP = {'a', 'an', 'the', 'on', 'of', 'in', 'and', 'for', 'to', 'with', 'by', 'at', 'from', 'as', 'is', 'are',
        'its', 'via', 'into', 'new', 'using', 'toward', 'towards'}


def make_key(msg):
    fam = 'anon'
    for a in msg.get('author', []) or []:
        if a.get('family'):
            fam = a['family']; break
        if a.get('name'):
            fam = a['name'].split()[0]; break
    fam = re.sub(r'[^a-z]', '', _ascii(fam).lower()) or 'anon'
    year = _year(msg) or 0
    title = _ascii(re.sub(r'<[^>]+>', '', (msg.get('title') or [''])[0])).lower()
    words = [w for w in re.findall(r'[a-z0-9]+', title) if w not in STOP]
    return f"{fam}{year}{words[0] if words else 'x'}"


def _year(msg):
    for k in ('published-print', 'published-online', 'issued', 'created'):
        dp = (msg.get(k) or {}).get('date-parts')
        if dp and dp[0] and dp[0][0]:
            return int(dp[0][0])
    return None


def summary(msg):
    auth = msg.get('author', []) or []
    names = ', '.join(f"{a.get('family', a.get('name', '?'))} {''.join(x[0] for x in (a.get('given') or '').split() if x)}"
                      for a in auth[:3]) + (' et al.' if len(auth) > 3 else '')
    return {
        'doi': norm_doi(msg['DOI']),
        'key': make_key(msg),
        'title': re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', '', (msg.get('title') or [''])[0])).strip(),
        'authors': names,
        'year': _year(msg),
        'journal': (msg.get('container-title') or [''])[0],
        'volume': msg.get('volume', ''),
        'issue': msg.get('issue', ''),
        'pages': msg.get('page', '') or msg.get('article-number', ''),
        'type': msg.get('type', ''),
        'publisher': msg.get('publisher', ''),
    }


def integrity(doi, msg, refresh=False):
    """-> (refusals, concerns, updates): retraction-type flags, expressions of concern, other updates (info)."""
    doi = norm_doi(doi)
    refusals, concerns, updates = [], [], []

    def classify(typ, where, notice):
        t = (typ or '').lower().replace('-', '_').replace(' ', '_')
        if not t:
            return
        tag = f'{where}: {t}' + (f' ({notice})' if notice and norm_doi(notice) != doi else '')
        if t in REFUSE_TYPES or 'retract' in t:
            refusals.append(tag)
        elif t in CONCERN_TYPES or 'concern' in t:
            concerns.append(tag)
        else:
            updates.append({'type': t, 'doi': norm_doi(notice) if notice else '', 'via': where})

    oa = openalex(doi, refresh)
    if oa and oa.get('is_retracted'):
        refusals.append('OpenAlex: is_retracted')
    for u in msg.get('update-to', []) or []:
        # A notice record: update-to names what it updates. It flags this DOI only if it points at itself
        # (Elsevier in-place retractions do exactly that).
        if norm_doi(u.get('DOI', '')) == doi:
            classify(u.get('type'), 'Crossref update-to', u.get('DOI'))
    for u in msg.get('updated-by', []) or []:
        classify(u.get('type'), f"Crossref updated-by[{u.get('source', '?')}]", u.get('DOI'))
    for item in update_notices(doi, refresh):
        for u in item.get('update-to', []) or []:
            if norm_doi(u.get('DOI', '')) == doi:
                classify(u.get('type'), 'Crossref notice', item.get('DOI'))
    for rel_type in (msg.get('relation') or {}):
        if 'retract' in rel_type.lower():
            refusals.append(f'Crossref relation: {rel_type}')
    title = re.sub(r'<[^>]+>', '', ((msg.get('title') or [''])[0] or '')).strip().lower()
    if title.startswith(('retracted', 'retraction', 'withdrawn')):
        refusals.append('title says retracted/withdrawn')
    # de-duplicate, keep order
    uniq = lambda xs: list(dict.fromkeys(xs))
    seen, upd = set(), []
    for u in updates:
        k = (u['type'], u['doi'])
        if k not in seen:
            seen.add(k); upd.append(u)
    return uniq(refusals), uniq(concerns), upd


def retraction_flags(doi, msg):
    """Backwards-compatible helper (cqc-dft name): retraction-type flags only."""
    return integrity(doi, msg)[0]


def _search_openalex(a):
    params = {'search': a.query, 'per-page': str(min(a.n, 50)),
              'filter': f'from_publication_date:{a.date_from}-01-01,to_publication_date:{a.date_to}-12-31,has_doi:true'}
    if a.sort:
        params['sort'] = a.sort
    d = json.loads(_get('https://api.openalex.org/works?' + urllib.parse.urlencode(params),
                        min_interval=0.15, tries=2, no_wait_429=True))
    out = []
    for w in d.get('results', []):
        out.append(dict(doi=norm_doi(w.get('doi') or ''), year=w.get('publication_year'),
                        venue=((w.get('primary_location') or {}).get('source') or {}).get('display_name') or '',
                        cites=w.get('cited_by_count'), title=w.get('title') or '', abstract=abstract_text(w),
                        retracted=bool(w.get('is_retracted'))))
    return out, d.get('meta', {}).get('count')


def _search_crossref(a):
    params = {'query.bibliographic': a.query, 'rows': str(min(a.n, 50)),
              'filter': f'from-pub-date:{a.date_from}-01-01,until-pub-date:{a.date_to}-12-31,type:journal-article',
              'select': 'DOI,title,container-title,issued,is-referenced-by-count,abstract'}
    if a.sort and a.sort.startswith('cited_by_count'):
        params.update(sort='is-referenced-by-count', order='desc')
    d = json.loads(_get('https://api.crossref.org/works?' + urllib.parse.urlencode(params), tries=4))['message']
    out = []
    for w in d.get('items', []):
        ab = re.sub(r'<[^>]+>', ' ', w.get('abstract') or '')
        out.append(dict(doi=norm_doi(w['DOI']), year=_year(w), venue=(w.get('container-title') or [''])[0],
                        cites=w.get('is-referenced-by-count'), title=re.sub(r'<[^>]+>', '', (w.get('title') or [''])[0]),
                        abstract=re.sub(r'\s+', ' ', ab).strip(), retracted=False))
    return out, d.get('total-results')


def _search_s2(a):
    params = {'query': a.query, 'limit': str(min(a.n, 50)), 'year': f'{a.date_from}-{a.date_to}',
              'fields': 'title,year,venue,citationCount,externalIds,abstract'}
    d = json.loads(_get('https://api.semanticscholar.org/graph/v1/paper/search?' + urllib.parse.urlencode(params),
                        min_interval=1.1, tries=2, no_wait_429=True))
    out = []
    for w in d.get('data', []) or []:
        doi = (w.get('externalIds') or {}).get('DOI')
        if not doi:
            continue
        out.append(dict(doi=norm_doi(doi), year=w.get('year'), venue=w.get('venue') or '', cites=w.get('citationCount'),
                        title=w.get('title') or '', abstract=w.get('abstract') or '', retracted=False))
    return out, d.get('total')


def cmd_search(a):
    order = {'auto': ['openalex', 'crossref', 's2'], 'openalex': ['openalex'], 'crossref': ['crossref'], 's2': ['s2']}[a.backend]
    fn = {'openalex': _search_openalex, 'crossref': _search_crossref, 's2': _search_s2}
    last = None
    for b in order:
        try:
            res, total = fn[b](a)
            break
        except (RateLimited, RuntimeError) as e:
            last = e
            print(f'-- {b} unavailable ({type(e).__name__}); trying next backend', file=sys.stderr)
    else:
        sys.exit(f'all search backends failed: {last}')
    for w in res:
        flag = ' [RETRACTED]' if w['retracted'] else ''
        print(f"{w['doi']} | {w['year']} | {w['venue'][:45]} | cites {w['cites']} | {w['title'][:150]}{flag}")
        if a.abstracts and w['abstract']:
            print('    ' + w['abstract'][:700])
    print(f'-- backend {b}; {total} hits in total')


def _fmt_updates(upd):
    return ', '.join(f"{u['type']}{' ' + u['doi'] if u['doi'] else ''}" for u in upd)


def cmd_check(a):
    bad = 0
    for doi in a.dois:
        msg = crossref(doi)
        if not msg:
            print(f'NOT FOUND in Crossref: {doi}')
            bad += 1
            continue
        s = summary(msg)
        ref, conc, upd = integrity(doi, msg)
        bad += bool(ref)
        print(f"{'FLAGGED' if ref else 'OK'} {s['doi']} [{s['key']}] {s['year']} {s['journal']} {s['volume']}"
              f"{'(' + s['issue'] + ')' if s['issue'] else ''} {s['pages']} | {s['authors']} | {s['title']}"
              + (f"  !! {'; '.join(ref)}" if ref else '') + (f"  ! concern: {'; '.join(conc)}" if conc else '')
              + (f"  (updates: {_fmt_updates(upd)})" if upd else ''))
    sys.exit(1 if bad else 0)


def cmd_abstract(a):
    oa = openalex(a.doi)
    if not oa:
        print('not in OpenAlex'); sys.exit(1)
    print(oa.get('title'))
    print(abstract_text(oa) or '(no abstract available in OpenAlex; read the landing page with WebFetch if needed)')


def _load():
    if not CAND.exists():
        return []
    return [json.loads(l) for l in CAND.read_text().splitlines() if l.strip()]


def _save(rows):
    CAND.parent.mkdir(parents=True, exist_ok=True)
    tmp = CAND.with_suffix('.tmp')
    tmp.write_text(''.join(json.dumps(r, ensure_ascii=False) + '\n' for r in rows))
    tmp.replace(CAND)


def _locked():
    p = CAND.with_suffix('.lock')
    p.parent.mkdir(parents=True, exist_ok=True)
    fh = open(p, 'a+')
    fcntl.flock(fh, fcntl.LOCK_EX)
    return fh


def _today():
    return datetime.date.today().isoformat()


def add_one(doi, topic, claim, relevance=3, section='', must_cite=False, force=False, allow_concern=False):
    """Verify one DOI and add it (or a further claim) to lit/candidates.jsonl. -> (ok, message)"""
    msg = crossref(doi)
    if not msg:
        return False, f'REFUSED: {doi} does not resolve in Crossref; it cannot be cited'
    s = summary(msg)
    ref, conc, upd = integrity(doi, msg)
    if ref:
        return False, f"REFUSED: {s['doi']} is flagged ({'; '.join(ref)})"
    if conc and not allow_concern:
        return False, f"REFUSED: {s['doi']} has an expression of concern ({'; '.join(conc)}); --allow-concern to override"
    if not csl(doi):
        return False, f"REFUSED: no CSL record for {s['doi']}"
    with _locked():
        rows = _load()
        for r in rows:
            if r['doi'] == s['doi']:
                r.setdefault('claims', [])
                if claim not in r['claims']:
                    r['claims'].append(claim)
                if r.get('status') == 'rejected' and force:
                    r['status'] = 'candidate'
                if must_cite:
                    r['must_cite'] = True
                r['updates'] = upd
                r['verified_on'] = _today()
                _save(rows)
                return True, f"already present as [{r['key']}]; claim added"
        keys = {r['key'] for r in rows}
        key = s['key']
        if key in keys:
            n = 2
            while f'{key}{chr(96 + n)}' in keys:
                n += 1
            key = f'{key}{chr(96 + n)}'
        s.update(key=key, topic=topic, claims=[claim], relevance=relevance, section=section,
                 status='candidate', must_cite=bool(must_cite), updates=upd, concerns=conc,
                 abstract_checked=bool(abstract_text(openalex(doi) or {})), verified_on=_today())
        rows.append(s)
        _save(rows)
    note = f"  (updates: {_fmt_updates(upd)})" if upd else ''
    return True, f"added [{key}] {s['year']} {s['journal']} | {s['title']}{note}"


def cmd_add(a):
    ok, m = add_one(a.doi, a.topic, a.claim, a.relevance, a.section, a.must_cite, a.force, a.allow_concern)
    print(m)
    if not ok:
        sys.exit(1)


def cmd_add_many(a):
    items = json.loads(Path(a.file).read_text())
    n_ok = n_bad = 0
    for it in items:
        ok, m = add_one(it['doi'], it['topic'], it['claim'], it.get('relevance', 3), it.get('section', ''),
                        it.get('must_cite', False), a.force, a.allow_concern)
        print(('  ' if ok else '!! ') + f"{norm_doi(it['doi'])}: {m}", flush=True)
        n_ok += ok; n_bad += not ok
    print(f'-- {n_ok} added/updated, {n_bad} refused')
    sys.exit(1 if n_bad else 0)


def cmd_drop(a):
    doi = norm_doi(a.doi)
    with _locked():
        rows = _load()
        hit = [r for r in rows if r['doi'] == doi]
        if not hit:
            sys.exit('not a candidate')
        for r in hit:
            r['status'] = 'rejected'; r['reject_reason'] = a.reason
        _save(rows)
    print(f'rejected {doi}')


def cmd_claim(a):
    """Replace one claim (exact text) of the candidate with key a.key; nothing else of the record changes."""
    with _locked():
        rows = _load()
        hit = [r for r in rows if r['key'] == a.key]
        if len(hit) != 1:
            sys.exit(f'{a.key}: {len(hit)} candidates with this key')
        r = hit[0]
        claims = r.get('claims', [])
        if a.old not in claims:
            sys.exit(f'{a.key}: no claim with this exact text; its claims are:\n  ' + '\n  '.join(claims))
        r['claims'] = [a.new if c == a.old else c for c in claims]
        _save(rows)
    print(f'[{a.key}] claim replaced')


def cmd_list(a):
    for r in sorted(_load(), key=lambda r: (r.get('topic', ''), -(r.get('year') or 0))):
        if r.get('status') == 'rejected' and not a.all:
            continue
        mc = ' MUST-CITE' if r.get('must_cite') else ''
        print(f"[{r['key']}]{mc} {r['year']} {r['journal'][:40]} | {r.get('topic')} | rel {r.get('relevance')} | "
              f"{r['title'][:90]} || {' / '.join(r.get('claims', []))[:200]}")


def _cited_keys(path, printed_only=False):
    """Citation keys in a manuscript. printed_only=True drops HTML comments first: the skeleton's planning notes
    (<!-- ... [@key] ... -->) are never printed, so a must-cite key that survives only in a note is NOT cited
    (review fix, 28 Sep 2026). The bibliography itself still scans everything, which is harmless (pandoc prints
    only the keys cited in the text) and keeps unknown keys in notes from going unnoticed."""
    p = Path(path)
    if not p.exists():
        return set()
    t = p.read_text()
    if printed_only:
        t = re.sub(r'<!--.*?-->', ' ', t, flags=re.S)
    # an "@" after a letter or digit is an e-mail address of the author block, not a citation (pandoc reads it so too)
    return set(re.findall(r'(?<![A-Za-z0-9])@([A-Za-z0-9_:-]+)', t))


def cmd_stats(a):
    from collections import Counter
    rows = [r for r in _load() if r.get('status') != 'rejected']
    print(f'{len(rows)} active candidates')
    yc = Counter(r['year'] for r in rows)
    recent = sum(v for k, v in yc.items() if k and k >= 2020)
    print(f'2020 or later: {recent} ({100 * recent / max(1, len(rows)):.0f}%)')
    print('by year:', dict(sorted(yc.items(), key=lambda kv: kv[0] or 0)))
    print('by topic:', dict(Counter(r.get('topic') for r in rows).most_common()))
    print('by section:', dict(Counter(r.get('section') or '-' for r in rows).most_common()))
    print('top journals:', Counter(r['journal'] for r in rows).most_common(12))
    cited = _cited_keys(MAIN_MS, printed_only=True)
    for r in rows:
        if r.get('must_cite'):
            print(f"must-cite [{r['key']}] {r['doi']}: {'cited' if r['key'] in cited else 'NOT YET CITED'} "
                  f"in {MAIN_MS.relative_to(REPO)}")
    with_upd = [r for r in rows if r.get('updates')]
    for r in with_upd:
        print(f"has published updates [{r['key']}]: {_fmt_updates(r['updates'])}")


def cmd_verify(a):
    """Fresh Crossref + OpenAlex + notice lookup for every active candidate."""
    rows = _load()
    bad = 0
    for r in rows:
        if r.get('status') == 'rejected':
            continue
        msg = crossref(r['doi'], refresh=not a.no_refresh)
        if not msg:
            print(f"!! [{r['key']}] {r['doi']}: no longer resolves in Crossref"); bad += 1; continue
        ref, conc, upd = integrity(r['doi'], msg, refresh=not a.no_refresh)
        r['updates'], r['concerns'], r['verified_on'] = upd, conc, _today()
        if ref:
            r['flags'] = ref
            print(f"!! [{r['key']}] {r['doi']}: FLAGGED {'; '.join(ref)}"); bad += 1
        elif conc:
            print(f"!  [{r['key']}] {r['doi']}: expression of concern {'; '.join(conc)}")
        else:
            print(f"ok [{r['key']}] {r['doi']}" + (f"  (updates: {_fmt_updates(upd)})" if upd else ''), flush=True)
    with _locked():
        cur = {x['doi']: x for x in _load()}
        for r in rows:
            if r['doi'] in cur:
                cur[r['doi']].update({k: r[k] for k in ('updates', 'concerns', 'verified_on', 'flags') if k in r})
        _save(list(cur.values()))
    n = sum(1 for r in rows if r.get('status') != 'rejected')
    print(f'-- verified {n} active candidates on {_today()}: {bad} problem(s)')
    sys.exit(1 if bad else 0)


def _cite_line(r):
    iss = f"({r['issue']})" if r.get('issue') else ''
    return (f"{r['authors']} ({r['year']}). {r['title']}. *{r['journal']}* {r.get('volume', '')}{iss}"
            f"{', ' + r['pages'] if r.get('pages') else ''}. doi:{r['doi']}")


def cmd_table(a):
    rows = [r for r in _load() if r.get('status') != 'rejected']
    recent = sum(1 for r in rows if (r.get('year') or 0) >= 2020)
    out = ['# SAPH literature pool (generated by `paper/_shared/tools/lit.py table`; do not edit by hand)', '',
           f'{len(rows)} active candidates, {recent} from 2020 onwards ({100 * recent / max(1, len(rows)):.0f} %). '
           'Every entry resolved in Crossref and passed the retraction check when it was added '
           '(`verified_on`); `lit.py verify` repeats the check. Cite as `[@key]`.', '']
    mc = [r for r in rows if r.get('must_cite')]
    if mc:
        out += ['## Must cite', ''] + [f"- `{r['key']}`: {_cite_line(r)}" for r in mc] + ['']
    by_topic = {}
    for r in rows:
        by_topic.setdefault(r.get('topic') or '-', []).append(r)
    for t in sorted(by_topic):
        out += [f'## {t}', '']
        for r in sorted(by_topic[t], key=lambda r: (-(r.get('year') or 0), r['key'])):
            upd = f" Published updates: {_fmt_updates(r['updates'])}." if r.get('updates') else ''
            out.append(f"- `{r['key']}`{' (must cite)' if r.get('must_cite') else ''}: {_cite_line(r)}. "
                       f"Use: {' / '.join(r.get('claims', []))}{upd}")
        out.append('')
    (LIT / 'REFS.md').write_text('\n'.join(out))
    print(f"wrote {LIT / 'REFS.md'} ({len(rows)} entries)")


def cmd_bib(a):
    """CSL-JSON for everything the manuscript cites."""
    rows = [r for r in _load() if r.get('status') != 'rejected']
    ms = Path(a.manuscript) if a.manuscript else MAIN_MS
    cited = _cited_keys(ms) if ms.exists() else None
    out, missing, seen = [], [], set()
    # Metadata corrections for records Crossref gets wrong (missing authors, mangled capitalisation):
    #   lit/csl_overrides.json  {"key": {"author": [{"family": "...", "given": "..."}], "title": "...", "_why": "..."}}
    # Every override must carry a "_why" saying where the correct value was read.
    overrides = {}
    src = LIT / 'csl_overrides.json'
    if src.exists():
        for k, v in json.loads(src.read_text()).items():
            if k.startswith('_'):
                continue
            if not v.get('_why'):
                sys.exit(f'{src}: override for {k} has no "_why" (evidence) field')
            overrides.setdefault(k, {}).update(v)
    for r in rows:
        if r['key'] in seen or (cited is not None and r['key'] not in cited):
            continue
        seen.add(r['key'])
        c = csl(r['doi'])
        if not c:
            missing.append(r['doi']); continue
        c = dict(c); c['id'] = r['key']
        c.pop('reference', None)          # Crossref's own reference list: not needed, and it bloats the file
        for f, v in overrides.get(r['key'], {}).items():
            if not f.startswith('_'):
                c[f] = v
        out.append(c)
    if cited is not None:
        unknown = sorted(k for k in cited if k not in seen and not k.startswith(('fig:', 'tbl:', 'eq:', 'sec:')))
        if unknown:
            print('cited in manuscript but not a verified candidate:', ', '.join(unknown))
        if ms.resolve() == MAIN_MS.resolve():
            printed = _cited_keys(ms, printed_only=True)
            for r in rows:
                if r.get('must_cite') and r['key'] not in printed:
                    print(f"WARNING: must-cite reference [{r['key']}] ({r['doi']}) is not cited yet")
    dest = Path(a.out) if a.out else ms.parent / 'references.json'
    dest.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(f'wrote {len(out)} references to {dest}' + (f'; {len(missing)} without CSL' if missing else ''))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest='cmd', required=True)
    s = sp.add_parser('search'); s.add_argument('query'); s.add_argument('--from', dest='date_from', type=int, default=2020)
    s.add_argument('--to', dest='date_to', type=int, default=2026); s.add_argument('--n', type=int, default=25)
    s.add_argument('--backend', choices=['auto', 'openalex', 'crossref', 's2'], default='auto')
    s.add_argument('--abstracts', action='store_true'); s.add_argument('--sort', default=None,
                                                                         help='e.g. cited_by_count:desc, publication_date:desc')
    s.set_defaults(f=cmd_search)
    s = sp.add_parser('check'); s.add_argument('dois', nargs='+'); s.set_defaults(f=cmd_check)
    s = sp.add_parser('abstract'); s.add_argument('doi'); s.set_defaults(f=cmd_abstract)
    s = sp.add_parser('add'); s.add_argument('doi'); s.add_argument('--topic', required=True)
    s.add_argument('--claim', required=True); s.add_argument('--relevance', type=int, default=3)
    s.add_argument('--section', default=''); s.add_argument('--must-cite', action='store_true')
    s.add_argument('--force', action='store_true', help='re-activate a rejected candidate')
    s.add_argument('--allow-concern', action='store_true', help='accept a DOI with an expression of concern')
    s.set_defaults(f=cmd_add)
    s = sp.add_parser('add-many'); s.add_argument('file'); s.add_argument('--force', action='store_true')
    s.add_argument('--allow-concern', action='store_true'); s.set_defaults(f=cmd_add_many)
    s = sp.add_parser('drop'); s.add_argument('doi'); s.add_argument('--reason', required=True)
    s.set_defaults(f=cmd_drop)
    s = sp.add_parser('claim'); s.add_argument('key'); s.add_argument('--old', required=True)
    s.add_argument('--new', required=True); s.set_defaults(f=cmd_claim)
    s = sp.add_parser('list'); s.add_argument('--all', action='store_true'); s.set_defaults(f=cmd_list)
    s = sp.add_parser('stats'); s.set_defaults(f=cmd_stats)
    s = sp.add_parser('verify'); s.add_argument('--no-refresh', action='store_true',
                                                help='use cached records (offline consistency check only)')
    s.set_defaults(f=cmd_verify)
    s = sp.add_parser('table'); s.set_defaults(f=cmd_table)
    s = sp.add_parser('bib')
    s.add_argument('--manuscript', help='cite-scan this file (default paper/main/manuscript.md; the SI: paper/si/manuscript.md)')
    s.add_argument('--out', help='write here instead of <manuscript dir>/references.json'); s.set_defaults(f=cmd_bib)
    a = p.parse_args()
    a.f(a)


if __name__ == '__main__':
    main()
