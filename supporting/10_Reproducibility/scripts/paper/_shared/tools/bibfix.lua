--[[ bibfix.lua -- tidy the citeproc output; run AFTER --citeproc.

DOIs in the reference list. The Elsevier CSL prints the DOI as the literal prefix "https://doi.org/" followed by a
link whose text is the bare DOI, so only half of the URL is clickable (and coloured). Merge the two into one link that
shows and points to the full https://doi.org/... address. The ACS CSL prints the full address as one link already;
it gets the same LaTeX treatment (below) so that the PDF never breaks inside "https://doi.org/".

LaTeX: the link is written as \href{URL}{\mbox{https://doi.org/10.xxxx/}\nolinkurl{SUFFIX}}, so a line can
break inside the DOI suffix (after / . - etc., see \UrlBreaks in manuscript.latex) or before the address, but
never inside "https://doi.org/10.xxxx/" (a plain \url broke it as "https:" / "//doi.org/..."). The escaping is
pandoc's own: the pieces are cut out of what pandoc's LaTeX writer produces for \url and \href.

Superscript citations (ACS profile, added 29 Sep 2026 for ACS Omega). citeproc renders an ACS citation as a
Superscript ("salt^1-3^"), placed where the [@key] stood, i.e. before any punctuation that follows.
  * A full stop or comma that directly follows the citation is moved in front of it, as ACS prints it
    ("amino sugars.^5^", "... salts,^3,4^ and"); colons, semicolons and closing brackets stay after it.
  * DOCX only: the citation is wrapped in the character style "Citation Number" (pandoc custom-style), so that
    docx_check.py can tell citation numbers from other superscripts (cm⁻¹, R²₂(8), Mg²⁺) and check that they run
    1..N in the order of first appearance.
Bracketed (Elsevier) citations are left as they are.
]]

local PREFIX = 'https://doi.org/'
local IS_LATEX = FORMAT:match('latex') ~= nil
local IS_DOCX = FORMAT:match('docx') ~= nil
local CITE_STYLE = 'Citation Number'

local function latex_of(inline)
  return pandoc.write(pandoc.Pandoc { pandoc.Plain { inline } }, 'latex'):gsub('%s+$', '')
end

local function plain_escape(s)
  -- pattern-escape a literal string for string.find/sub
  return (s:gsub('[%^%$%(%)%%%.%[%]%*%+%-%?]', '%%%0'))
end

local function doi_link(target, title, attr, shown)
  -- shown: the DOI as citeproc printed it (not percent-encoded like the target, e.g. SICI DOIs with < >)
  local text = (shown and shown ~= '') and (PREFIX .. shown) or target
  local link = pandoc.Link({ pandoc.Str(text) }, target, title, attr)
  if not IS_LATEX then return link end
  -- \url{https://doi.org/<doi, escaped for \url>}
  local u = latex_of(link)
  local url_doi = u:match('^\\url{' .. plain_escape(PREFIX) .. '(.*)}$')
  -- \href{<url, escaped for \href>}{x}
  local h = latex_of(pandoc.Link({ pandoc.Str('x') }, target, title))
  local href_url = h:match('^\\href{(.*)}{x}$')
  if not url_doi or not href_url then return link end        -- unexpected writer output: keep pandoc's
  -- the registrant ("10.1016/") stays with the prefix; the line may break only inside the suffix
  url_doi = url_doi:gsub('\\%%3[Cc]', '<'):gsub('\\%%3[Ee]', '>')   -- show SICI brackets, link keeps %3C/%3E
  local reg, suffix = url_doi:match('^(10%.[%d%.]+/)(.+)$')
  if not reg then reg, suffix = '', url_doi end
  return pandoc.RawInline('latex',
    '\\href{' .. href_url .. '}{\\mbox{' .. PREFIX .. reg .. '}\\nolinkurl{' .. suffix .. '}}')
end

local function merge(inlines)
  local out = pandoc.Inlines {}
  local i = 1
  local changed = false
  while i <= #inlines do
    local a, b = inlines[i], inlines[i + 1]
    if a.t == 'Str' and b and b.t == 'Link' and a.text:sub(-#PREFIX) == PREFIX
        and b.target:sub(1, #PREFIX):lower() == PREFIX then
      local before = a.text:sub(1, #a.text - #PREFIX)
      if before ~= '' then out:insert(pandoc.Str(before)) end
      out:insert(doi_link(b.target, b.title, b.attr, pandoc.utils.stringify(b.content)))
      i = i + 2
      changed = true
    elseif IS_LATEX and a.t == 'Link' and a.target:sub(1, #PREFIX):lower() == PREFIX
        and pandoc.utils.stringify(a.content):sub(1, #PREFIX):lower() == PREFIX then
      -- ACS: the whole address is one link whose text is the address
      out:insert(doi_link(a.target, a.title, a.attr, pandoc.utils.stringify(a.content):sub(#PREFIX + 1)))
      i = i + 1
      changed = true
    else
      out:insert(a)
      i = i + 1
    end
  end
  return changed and out or nil
end

local function is_sup_cite(el)
  return el.t == 'Cite' and #el.content > 0 and el.content[1].t == 'Superscript'
end

-- superscript citations: punctuation in front, and the DOCX character style
local function sup_cites(inlines)
  local out = pandoc.Inlines {}
  local changed = false
  local i = 1
  while i <= #inlines do
    local el = inlines[i]
    if is_sup_cite(el) then
      local nxt = inlines[i + 1]
      if nxt and nxt.t == 'Str' then
        local p = nxt.text:sub(1, 1)
        if p == '.' or p == ',' then
          out:insert(pandoc.Str(p))
          local rest = nxt.text:sub(2)
          if rest ~= '' then inlines[i + 1] = pandoc.Str(rest) else inlines:remove(i + 1) end
          changed = true
        end
      end
      if IS_DOCX then
        el.content = pandoc.Inlines { pandoc.Span(el.content, { ['custom-style'] = CITE_STYLE }) }
        changed = true
      end
      out:insert(el)
    else
      out:insert(el)
    end
    i = i + 1
  end
  return changed and out or nil
end

function Pandoc(doc)
  doc.blocks = doc.blocks:walk {
    Div = function(d)
      if d.identifier == 'refs' then
        return d:walk { Inlines = merge }, false
      end
    end,
  }
  doc.blocks = doc.blocks:walk { Inlines = sup_cites }
  return doc
end
