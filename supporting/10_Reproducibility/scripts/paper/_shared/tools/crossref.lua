--[[ crossref.lua -- figure, table, equation and section cross-references for the manuscripts.

Run BEFORE --citeproc:  pandoc ... --lua-filter crossref.lua --citeproc ...

Labels (numbered in order of appearance, separately for each kind):
  ![Caption.](figures/x.png){#fig:label}          figure
  Table: Caption. {#tbl:label}                    table (caption line of a pipe table)
  $$ ... $$ {#eq:label}                           display equation ("$$...$$ {#eq:x}," keeps the comma)
  # 2.1. Heading {#sec:label}                     section; the number is read from the heading text

References (plain text in DOCX, hyperlinks in PDF):
  @fig:a  or [@fig:a]           -> Fig. 1          [@fig:a; @fig:b]      -> Figs. 1 and 2
  [@fig:a; @fig:b; @fig:c]      -> Figs. 1-3       [@fig:a a; @fig:a b]  -> Figs. 2a and 2b (panel suffix)
  @tbl:x -> Table 2     @eq:x -> Eq. (3)     [@eq:x; @eq:y] -> Eqs. (3) and (4)     @sec:x -> Section 2.1
  -@fig:a -> 1 (number only)   [see @fig:a] -> see Fig. 1
  @fig:a [@key1; @key2] -> Fig. 1 [4,5]   (pandoc parses this as one Cite; it is split, and citeproc numbers
                                            the bibliography part)
Captions become "Fig. 1. ..." and "Table 1. ..." (DOCX: text inserted here; LaTeX: \setcounter before the float
so that LaTeX's own caption label carries exactly the same number).

Also:
  * Unicode super/subscript characters in text (cm⁻¹, H₂O) become real Superscript/Subscript elements, so they
    render in any font (TeX Gyre Termes has no U+207B) and look like Word's own superscripts.
  * Images without an explicit width get one: natural size from the PNG's pHYs chunk (300 dpi assumed when
    absent), capped at max-figure-width (default 16 cm) and max-figure-height (default 20 cm).
    pdf-max-figure-height (cm) lowers the height cap for LaTeX/PDF output only, so that a tall figure and its
    long caption fit on one PDF page; the Word file keeps the full size (final build, 29 Sep 2026).
    An image attribute pdf-max-height=<cm> lowers the cap further for that one figure, again for LaTeX/PDF only.
  * Front matter from the YAML block (author-block, abstract, keywords) is turned into body blocks, so Word and
    PDF show the same thing in the same order. author-block (Markdown blocks) is the title-page author block,
    printed in the "Author" style; the older one-line authors-placeholder still works. Neither is invented when
    the YAML has none (ACS adaptation, 29 Sep 2026). abstract-title sets the abstract heading (default Abstract).
  * An empty Div {#si-contents} is filled with a contents list: one line per numbered level-1 section
    ("S4. Harmonic frequencies ...: Tables S5 and S6; Fig. S1"), from the headings and the figures and tables
    that follow each (Supporting Information first page; DOCX style "Contents List").
  * Metadata crossref-map=<path> writes the label -> number map as JSON (used by the build checks).
  * Metadata xref-style=acs (build.py, profile acs-omega) switches to the labels of the ACS journals: "Figure 1" and
    "Figures 1 and 2" in the text and "Figure 1." in the caption, and "eq 1" / "eqs 1 and 2" without parentheses in
    the text (the equation itself keeps its "(1)"). Without it (profile elsevier) the labels are "Fig. 1", "Fig. 1."
    and "Eq. (1)" as described above (ACS integration review, 30 Sep 2026).
Any undefined or duplicate label, a mixed citation/cross-reference bracket, or an unconsumed {#...} label is a
hard error.
]]

local stringify = pandoc.utils.stringify
local IS_LATEX = FORMAT:match('latex') ~= nil or FORMAT:match('beamer') ~= nil
local IS_DOCX = FORMAT:match('docx') ~= nil
local NBSP = '\u{00A0}'

local KINDS = { fig = true, tbl = true, eq = true, sec = true }
local NAMES = {
  fig = { 'Fig.', 'Figs.' },
  tbl = { 'Table', 'Tables' },
  eq = { 'Eq.', 'Eqs.' },
  sec = { 'Section', 'Sections' },
}
local LATEX_COUNTER = { fig = 'figure', tbl = 'table' }
local FIG_LABEL = 'Fig.'   -- caption label; 'Figure' with xref-style=acs
local EQ_PAREN = true      -- "Eq. (1)"; false with xref-style=acs ("eq 1")

local labels = {}          -- label -> { kind=, num= (string), order= }
local counters = { fig = 0, tbl = 0, eq = 0 }
local order = {}
local errors = {}
local warnings = {}
local max_w_cm, max_h_cm = 16.0, 20.0
local PREFIX = ''          -- metadata number-prefix, e.g. 'S' for Supporting Information

local function err(msg) errors[#errors + 1] = msg end
local function warn(msg) warnings[#warnings + 1] = msg end

local function kind_of(label)
  local k = label:match('^(%a+):')
  if k and KINDS[k] then return k end
  return nil
end

local function define(label, kind, num)
  if labels[label] then
    err('duplicate label #' .. label)
    return
  end
  labels[label] = { kind = kind, num = tostring(num) }
  order[#order + 1] = label
end

---------------------------------------------------------------------------------------------------
-- Unicode superscripts / subscripts -> Superscript / Subscript elements
---------------------------------------------------------------------------------------------------
local SUP = {
  [0x2070] = '0', [0x00B9] = '1', [0x00B2] = '2', [0x00B3] = '3', [0x2074] = '4', [0x2075] = '5',
  [0x2076] = '6', [0x2077] = '7', [0x2078] = '8', [0x2079] = '9', [0x207A] = '+', [0x207B] = '\u{2212}',
  [0x207C] = '=', [0x207D] = '(', [0x207E] = ')', [0x2071] = 'i', [0x207F] = 'n',
}
local SUB = {
  [0x2080] = '0', [0x2081] = '1', [0x2082] = '2', [0x2083] = '3', [0x2084] = '4', [0x2085] = '5',
  [0x2086] = '6', [0x2087] = '7', [0x2088] = '8', [0x2089] = '9', [0x208A] = '+', [0x208B] = '\u{2212}',
  [0x208C] = '=', [0x208D] = '(', [0x208E] = ')', [0x2090] = 'a', [0x2091] = 'e', [0x2092] = 'o',
  [0x2093] = 'x', [0x2095] = 'h', [0x2096] = 'k', [0x2097] = 'l', [0x2098] = 'm', [0x2099] = 'n',
  [0x209A] = 'p', [0x209B] = 's', [0x209C] = 't',
}

local function split_scripts(s)
  -- fast path: nothing to do
  local has = false
  for _, c in utf8.codes(s) do
    if SUP[c] or SUB[c] then has = true; break end
  end
  if not has then return nil end
  local out, buf, mode, run = {}, {}, nil, {}
  local function flush_text()
    if #buf > 0 then out[#out + 1] = pandoc.Str(table.concat(buf)); buf = {} end
  end
  local function flush_run()
    if mode then
      local txt = table.concat(run)
      out[#out + 1] = (mode == 'sup') and pandoc.Superscript { pandoc.Str(txt) } or pandoc.Subscript { pandoc.Str(txt) }
      run, mode = {}, nil
    end
  end
  for _, c in utf8.codes(s) do
    local m = SUP[c] and 'sup' or (SUB[c] and 'sub' or nil)
    if m then
      if mode ~= m then flush_run(); flush_text(); mode = m end
      run[#run + 1] = SUP[c] or SUB[c]
    else
      flush_run()
      buf[#buf + 1] = utf8.char(c)
    end
  end
  flush_run(); flush_text()
  return out
end

---------------------------------------------------------------------------------------------------
-- Figure sizes
---------------------------------------------------------------------------------------------------
local function png_dpi(contents)
  -- exact pixels-per-metre from the pHYs chunk (pandoc.image.size truncates 299.99 to 299)
  if contents:sub(1, 8) ~= '\137PNG\r\n\26\n' then return nil end
  local pos = 9
  while pos + 8 <= #contents do
    local len, typ = string.unpack('>I4c4', contents, pos)
    if typ == 'pHYs' then
      local ppux, ppuy, unit = string.unpack('>I4I4B', contents, pos + 8)
      if unit == 1 and ppux > 0 then return ppux * 0.0254, ppuy * 0.0254 end
      return nil
    elseif typ == 'IDAT' or typ == 'IEND' then
      return nil
    end
    pos = pos + 12 + len
  end
  return nil
end

local function size_image(img)
  -- per-figure PDF height cap (cm): {#fig:x pdf-max-height=12}; LaTeX only, removed from every output
  -- (29 Sep 2026: Fig. 2 gained two panels and its caption no longer fitted one PDF page under the global cap)
  local fig_h_cm = max_h_cm
  if img.attributes['pdf-max-height'] then
    if IS_LATEX then fig_h_cm = math.min(max_h_cm, tonumber(img.attributes['pdf-max-height']) or max_h_cm) end
    img.attributes['pdf-max-height'] = nil
  end
  if img.attributes.width or img.attributes.height then
    local w = img.attributes.width
    if w then
      local v, u = w:match('^([%d%.]+)%s*(%a*)$')
      v = tonumber(v)
      local cm = v and ((u == 'cm' and v) or (u == 'mm' and v / 10) or (u == 'in' and v * 2.54)) or nil
      if cm and cm > max_w_cm + 1e-6 then
        warn(('%s: explicit width %s is wider than %.1f cm (kept)'):format(img.src, w, max_w_cm))
      end
    end
    return img
  end
  local ok, mt, contents = pcall(pandoc.mediabag.fetch, img.src)
  if not ok or not contents then
    err('cannot read image ' .. img.src)
    return img
  end
  local sz = pandoc.image.size(contents)
  local dx, dy = png_dpi(contents)
  if not dx then
    dx, dy = sz.dpi_horz, sz.dpi_vert
    if not dx or dx < 150 then
      warn(('%s has no usable resolution metadata (%s dpi); assuming 300 dpi'):format(img.src, tostring(dx)))
      dx, dy = 300, 300
    end
  end
  local w = sz.width / dx * 2.54
  local h = sz.height / dy * 2.54
  local scale = math.min(1.0, max_w_cm / w, fig_h_cm / h)
  img.attributes.width = ('%.2fcm'):format(w * scale)
  return img
end

---------------------------------------------------------------------------------------------------
-- Pass 1: collect labels in document order
---------------------------------------------------------------------------------------------------
local function eq_label_from(el)
  if el and el.t == 'Str' then
    local label, rest = el.text:match('^{#(eq:[^}%s]+)}(.*)$')
    if label then return label, rest end
    if el.text:match('^{#') then err('malformed equation label ' .. el.text) end
  end
  return nil
end

-- For a list of inlines, returns positions of (display math, label) pairs.
local function find_equations(inlines)
  local found = {}
  local i = 1
  while i <= #inlines do
    local el = inlines[i]
    if el.t == 'Math' and el.mathtype == 'DisplayMath' then
      local j = i + 1
      while inlines[j] and (inlines[j].t == 'Space' or inlines[j].t == 'SoftBreak') do j = j + 1 end
      local label, rest = eq_label_from(inlines[j])
      found[#found + 1] = { math = i, label_pos = label and j or nil, label = label, rest = rest }
      i = label and j + 1 or i + 1
    else
      i = i + 1
    end
  end
  return found
end

local function header_number(h)
  local first = h.content[1]
  if first and first.t == 'Str' then
    local n = first.text:match('^(%d[%d%.]*)$')
    if n then return (n:gsub('%.$', '')) end
  end
  return nil
end

local function collect(blocks)
  blocks:walk {
    traverse = 'topdown',
    Figure = function(fig)
      counters.fig = counters.fig + 1
      if fig.identifier ~= '' then
        if kind_of(fig.identifier) ~= 'fig' then
          err('figure label must start with fig: (#' .. fig.identifier .. ')')
        end
        define(fig.identifier, 'fig', counters.fig)
      else
        warn('figure ' .. counters.fig .. ' has no {#fig:...} label')
      end
    end,
    Table = function(tbl)
      local cap = stringify(tbl.caption.long or {})
      if cap == '' then
        if tbl.identifier ~= '' then err('table #' .. tbl.identifier .. ' has a label but no caption') end
        return
      end
      counters.tbl = counters.tbl + 1
      if tbl.identifier ~= '' then
        if kind_of(tbl.identifier) ~= 'tbl' then
          err('table label must start with tbl: (#' .. tbl.identifier .. ')')
        end
        define(tbl.identifier, 'tbl', counters.tbl)
      else
        warn('table ' .. counters.tbl .. ' has no {#tbl:...} label')
      end
    end,
    Para = function(p)
      for _, e in ipairs(find_equations(p.content)) do
        if e.label then
          counters.eq = counters.eq + 1
          define(e.label, 'eq', counters.eq)
        end
      end
    end,
    Header = function(h)
      if h.identifier:match('^sec:') then
        local n = header_number(h)
        if n then define(h.identifier, 'sec', n)
        else labels[h.identifier] = { kind = 'sec', num = nil } end
      end
    end,
  }
end

---------------------------------------------------------------------------------------------------
-- Reference text
---------------------------------------------------------------------------------------------------
local function is_panel_suffix(s)
  -- "a", "b", "a,b", "a-c", "(a)" : glued to the number, no space
  return s:match('^%(?[a-z]%)?$') or s:match('^[a-z][,%-\u{2013}][a-z]$') ~= nil
end

local function fmt_num(kind, num)
  if kind ~= 'sec' and num ~= '??' then num = PREFIX .. num end
  if kind == 'eq' and EQ_PAREN then return '(' .. num .. ')' end
  return num
end

-- items: list of { kind=, num=, suffix= } of ONE kind, in citation order -> list of strings
local function number_list(kind, items)
  local plain = true
  for _, it in ipairs(items) do
    if it.suffix ~= '' or not tonumber(it.num) then plain = false end   -- '??' (undefined) or panel suffix
  end
  local parts = {}
  if plain and kind ~= 'sec' then
    -- sort numerically and compress runs of three or more into a range
    local nums = {}
    for _, it in ipairs(items) do nums[#nums + 1] = tonumber(it.num) end
    table.sort(nums)
    local uniq = {}
    for _, n in ipairs(nums) do if uniq[#uniq] ~= n then uniq[#uniq + 1] = n end end
    local i = 1
    while i <= #uniq do
      local j = i
      while uniq[j + 1] and uniq[j + 1] == uniq[j] + 1 do j = j + 1 end
      if j - i >= 2 then
        parts[#parts + 1] = fmt_num(kind, tostring(uniq[i])) .. '\u{2013}' .. fmt_num(kind, tostring(uniq[j]))
        i = j + 1
      else
        parts[#parts + 1] = fmt_num(kind, tostring(uniq[i]))
        i = i + 1
      end
    end
  else
    for _, it in ipairs(items) do parts[#parts + 1] = fmt_num(kind, it.num) .. it.suffix end
  end
  return parts
end

local function join_and(parts)
  if #parts == 1 then return parts[1] end
  return table.concat(parts, ', ', 1, #parts - 1) .. ' and ' .. parts[#parts]
end

local function text_to_inlines(s)
  -- keep NBSP inside tokens; ordinary spaces become Space
  local out = {}
  for word in s:gmatch('[^ ]+') do
    if #out > 0 then out[#out + 1] = pandoc.Space() end
    out[#out + 1] = pandoc.Str(word)
  end
  return out
end

local xref_cite  -- forward declaration

local function Cite(el)
  local n_x = 0
  for _, c in ipairs(el.citations) do if kind_of(c.id) then n_x = n_x + 1 end end
  if n_x == 0 then return nil end
  if n_x < #el.citations then
    -- "@fig:a [@ref1; @ref2]" is parsed by pandoc as ONE Cite (in-text reference plus bracketed "locator"):
    -- split it into the cross-reference text and an ordinary bibliography Cite that citeproc numbers.
    local first = el.citations[1]
    local rest = pandoc.List()
    for i = 2, #el.citations do rest:insert(el.citations[i]) end
    local rest_is_bib = true
    for _, c in ipairs(rest) do
      if kind_of(c.id) or c.mode ~= 'NormalCitation' then rest_is_bib = false end
    end
    if kind_of(first.id) and first.mode == 'AuthorInText' and rest_is_bib then
      local head = xref_cite(pandoc.Cite({ pandoc.Str('@' .. first.id) }, { first }))
      local keys = {}
      for _, c in ipairs(rest) do keys[#keys + 1] = '@' .. c.id end
      local out = pandoc.List(head)
      out:insert(pandoc.Space())
      out:insert(pandoc.Cite({ pandoc.Str('[' .. table.concat(keys, '; ') .. ']') }, rest))
      return out
    end
    err('mixed bibliography citation and cross-reference in one bracket: ' .. stringify(el.content) ..
      ' (write the cross-reference and the citation separately, e.g. "Fig. 1 [@key]" as "@fig:x [@key]")')
    return nil
  end
  return xref_cite(el)
end

xref_cite = function(el)
  -- group consecutive citations of the same kind
  local groups = {}
  local suppress = false
  local prefix, trailing = nil, nil
  for idx, c in ipairs(el.citations) do
    local kind = kind_of(c.id)
    local lab = labels[c.id]
    if not lab then
      err('undefined cross-reference @' .. c.id)
      lab = { kind = kind, num = '??' }
    elseif lab.num == nil then
      err('@' .. c.id .. ' refers to a heading without a leading number')
      lab = { kind = kind, num = '??' }
    end
    if c.mode == 'SuppressAuthor' then suppress = true end
    if idx == 1 and #c.prefix > 0 then prefix = c.prefix end
    local sfx = stringify(c.suffix):gsub('^%s+', '')
    local suffix = ''
    if sfx ~= '' then
      if is_panel_suffix(sfx) then
        suffix = sfx
      elseif idx == #el.citations then
        trailing = c.suffix
      else
        err('unsupported suffix "' .. sfx .. '" in cross-reference @' .. c.id)
      end
    end
    local g = groups[#groups]
    if not g or g.kind ~= kind then
      g = { kind = kind, items = {} }
      groups[#groups + 1] = g
    end
    g.items[#g.items + 1] = { kind = kind, num = lab.num, suffix = suffix, id = c.id }
  end
  local chunks = {}
  for _, g in ipairs(groups) do
    local parts = number_list(g.kind, g.items)
    local nums = join_and(parts)
    if suppress then
      chunks[#chunks + 1] = nums
    else
      local single = (#g.items == 1)
      chunks[#chunks + 1] = NAMES[g.kind][single and 1 or 2] .. NBSP .. nums
    end
  end
  local text = join_and(chunks)
  local inl = text_to_inlines(text)
  if IS_LATEX and #el.citations == 1 then
    inl = { pandoc.Link(inl, '#' .. el.citations[1].id) }
  end
  local out = pandoc.List()
  if prefix then out:extend(prefix); out:insert(pandoc.Space()) end
  out:extend(inl)
  if trailing then out:extend(trailing) end
  return out
end

---------------------------------------------------------------------------------------------------
-- Pass 2 transformations
---------------------------------------------------------------------------------------------------
local function caption_prefix(kind, num, long)
  -- long: list of blocks (caption.long); prepend a bold "Fig. 1." to the first Plain/Para
  local lbl = (kind == 'fig' and FIG_LABEL or 'Table') .. NBSP .. PREFIX .. num .. '.'
  local first = long[1]
  if first and (first.t == 'Plain' or first.t == 'Para') then
    first.content:insert(1, pandoc.Space())
    first.content:insert(1, pandoc.Strong { pandoc.Str(lbl) })
  else
    long:insert(1, pandoc.Plain { pandoc.Strong { pandoc.Str(lbl) } })
  end
  return long
end

local fig_seen, tbl_seen = 0, 0

local function docx_id(id)
  -- Word bookmark names may contain only letters, digits and underscores
  return (id:gsub('[^%w_]', '_'))
end

local function Figure(fig)
  fig_seen = fig_seen + 1
  fig.content = fig.content:walk { Image = size_image }
  if IS_LATEX then
    return { pandoc.RawBlock('latex', ('\\setcounter{figure}{%d}'):format(fig_seen - 1)), fig }
  end
  fig.caption.long = caption_prefix('fig', fig_seen, fig.caption.long)
  if IS_DOCX then fig.identifier = docx_id(fig.identifier) end
  return fig
end

local function Table(tbl)
  if stringify(tbl.caption.long or {}) == '' then return nil end
  tbl_seen = tbl_seen + 1
  if IS_LATEX then
    -- PDF only: a table of 10 or more columns prints at 8 pt (\densetable, manuscript.latex)
    local dense = #tbl.colspecs >= 10 and '\\densetable' or ''
    return { pandoc.RawBlock('latex', ('\\setcounter{table}{%d}%s'):format(tbl_seen - 1, dense)), tbl }
  end
  tbl.caption.long = caption_prefix('tbl', tbl_seen, tbl.caption.long)
  if IS_DOCX then tbl.identifier = docx_id(tbl.identifier) end
  return tbl
end

local function latex_equation(tex, label, num)
  if label then
    return ('\\begin{linenomath*}\\begin{equation}\\tag{%s}\\label{%s}\n%s\n\\end{equation}\\end{linenomath*}')
        :format(PREFIX .. num, label, tex)
  end
  return ('\\begin{linenomath*}\\begin{equation*}\n%s\n\\end{equation*}\\end{linenomath*}'):format(tex)
end

local function docx_equation(math_el, label, num)
  local math_para = pandoc.Div({ pandoc.Para { math_el } }, { ['custom-style'] = 'Equation' })
  local num_para = pandoc.Div({ pandoc.Para { pandoc.Str('(' .. PREFIX .. num .. ')') } }, { ['custom-style'] = 'EquationNumber' })
  local row = pandoc.Row { pandoc.Cell { pandoc.Plain {} }, pandoc.Cell { math_para }, pandoc.Cell { num_para } }
  local colspecs = { { 'AlignLeft', 0.08 }, { 'AlignCenter', 0.84 }, { 'AlignRight', 0.08 } }
  local body = { attr = pandoc.Attr(), body = { row }, head = {}, row_head_columns = 0 }
  local tbl = pandoc.Table(pandoc.Caption(), colspecs, pandoc.TableHead(), { body }, pandoc.TableFoot())
  return pandoc.Div({ tbl }, pandoc.Attr(docx_id(label), { 'numbered-equation' }))
end

-- Split a paragraph around its display equations.
local function equation_para(p)
  local eqs = find_equations(p.content)
  if #eqs == 0 then return nil end
  local inl = p.content
  local blocks = pandoc.List()
  local cur = pandoc.List()
  local skip_to = 0
  local eqi = 1
  local function trim(list)
    while #list > 0 and (list[1].t == 'Space' or list[1].t == 'SoftBreak') do list:remove(1) end
    while #list > 0 and (list[#list].t == 'Space' or list[#list].t == 'SoftBreak') do list:remove(#list) end
    return list
  end
  local i = 1
  while i <= #inl do
    local e = eqs[eqi]
    if e and i == e.math then
      local m = inl[i]
      local tex = m.text
      local rest = e.rest or ''
      local punct, after = rest:match('^([%.,;:])(.*)$')
      if punct then tex = tex .. punct; rest = after end
      local num = e.label and labels[e.label].num or nil
      if IS_LATEX then
        cur:insert(pandoc.RawInline('latex', latex_equation(tex, e.label, num)))
      else
        trim(cur)
        if #cur > 0 then blocks:insert(pandoc.Para(cur)) end
        cur = pandoc.List()
        if e.label and IS_DOCX then
          blocks:insert(docx_equation(pandoc.Math('DisplayMath', tex), e.label, num))
        elseif e.label then
          blocks:insert(pandoc.Div({ pandoc.Para { pandoc.Math('DisplayMath', tex), pandoc.Space(),
            pandoc.Str('(' .. num .. ')') } }, pandoc.Attr(e.label, { 'numbered-equation' })))
        else
          blocks:insert(pandoc.Para { pandoc.Math('DisplayMath', tex) })
        end
      end
      if rest ~= '' then cur:insert(pandoc.Str(rest)) end
      i = (e.label_pos or e.math) + 1
      eqi = eqi + 1
    else
      cur:insert(inl[i])
      i = i + 1
    end
  end
  if IS_LATEX then
    return pandoc.Para(cur)
  end
  trim(cur)
  if #cur > 0 then blocks:insert(pandoc.Para(cur)) end
  return blocks
end

-- LaTeX: a paragraph that begins with a display equation (the usual Markdown layout: text, blank line, $$..$$)
-- would start a new TeX paragraph with an empty line above the equation. Append it to the preceding
-- paragraph instead; text after the equation keeps its own paragraph if the author left a blank line.
local function starts_with_display(b)
  if b.t ~= 'Para' then return false end
  for _, el in ipairs(b.content) do
    if el.t ~= 'Space' and el.t ~= 'SoftBreak' then
      return el.t == 'Math' and el.mathtype == 'DisplayMath'
    end
  end
  return false
end

local function join_display_paras(blocks)
  local out = pandoc.Blocks {}
  local changed = false
  for _, b in ipairs(blocks) do
    local prev = out[#out]
    if prev and prev.t == 'Para' and starts_with_display(b) then
      prev.content:insert(pandoc.SoftBreak())
      prev.content:extend(b.content)
      changed = true
    else
      out:insert(b)
    end
  end
  return changed and out or nil
end

local function Str(s)
  if s.text:match('{#%a+:') then
    err('unconsumed label "' .. s.text .. '" (labels go after an image, a table caption or $$...$$)')
  end
  return split_scripts(s.text)
end

local function Header(h)
  if IS_DOCX and h.identifier:match('[^%w_]') then
    -- Word bookmark names may not contain ':' or '-'; nothing links to headings in the DOCX anyway
    h.identifier = docx_id(h.identifier)
  end
  return h
end

---------------------------------------------------------------------------------------------------
-- Front matter and bibliography placement
---------------------------------------------------------------------------------------------------
local function front_matter(doc)
  local meta = doc.meta
  local blocks = pandoc.List()
  local ab_block = meta['author-block']
  local placeholder = meta['authors-placeholder'] and stringify(meta['authors-placeholder']) or ''
  if not meta.author and ab_block then
    local bl
    local ty = pandoc.utils.type(ab_block)
    if ty == 'Blocks' then bl = ab_block
    elseif ty == 'Inlines' then bl = pandoc.Blocks { pandoc.Para(ab_block) }
    else bl = pandoc.Blocks { pandoc.Para(text_to_inlines(stringify(ab_block))) } end
    bl = bl:walk { Str = function(s) return split_scripts(s.text) end }
    blocks:insert(pandoc.Div(bl, pandoc.Attr('', { 'authors' }, { ['custom-style'] = 'Author' })))
    meta['author-block'] = nil
  elseif not meta.author and placeholder ~= '' then
    blocks:insert(pandoc.Div({ pandoc.Para(text_to_inlines(placeholder)) },
      pandoc.Attr('', { 'authors' }, { ['custom-style'] = 'Author' })))
  end
  if meta.abstract then
    local atitle = meta['abstract-title'] and stringify(meta['abstract-title']) or 'Abstract'
    blocks:insert(pandoc.Header(1, text_to_inlines(atitle), pandoc.Attr('abstract', { 'unnumbered' })))
    local ab = meta.abstract
    local ty = pandoc.utils.type(ab)
    if ty == 'Blocks' then
      blocks:extend(ab)
    elseif ty == 'Inlines' then
      blocks:insert(pandoc.Para(ab))
    else
      blocks:insert(pandoc.Para { pandoc.Str(stringify(ab)) })
    end
    meta.abstract = nil
  end
  if meta.keywords then
    local kws = {}
    local kw = meta.keywords
    if pandoc.utils.type(kw) == 'List' then
      for _, k in ipairs(kw) do kws[#kws + 1] = stringify(k) end
    else
      kws[1] = stringify(kw)
    end
    local inl = pandoc.List { pandoc.Strong { pandoc.Str('Keywords:') }, pandoc.Space() }
    inl:extend(text_to_inlines(table.concat(kws, '; ')))
    blocks:insert(pandoc.Div({ pandoc.Para(inl) }, pandoc.Attr('', { 'keywords' }, { ['custom-style'] = 'Keywords' })))
  end
  local has_refs = false
  doc.blocks:walk { Div = function(d) if d.identifier == 'refs' then has_refs = true end end }
  local body = blocks
  body:extend(doc.blocks)
  if not has_refs then
    body:insert(pandoc.Header(1, { pandoc.Str('References') }, pandoc.Attr('references', { 'unnumbered' })))
    body:insert(pandoc.Div({}, pandoc.Attr('refs')))
  end
  doc.blocks = pandoc.Blocks(body)
  return doc
end

---------------------------------------------------------------------------------------------------
-- Contents list (Supporting Information): Div {#si-contents}
---------------------------------------------------------------------------------------------------
local function float_list(kind, nums)
  if #nums == 0 then return nil end
  local items = {}
  for _, n in ipairs(nums) do items[#items + 1] = { kind = kind, num = tostring(n), suffix = '' } end
  return NAMES[kind][#nums == 1 and 1 or 2] .. NBSP .. join_and(number_list(kind, items))
end

local function contents_list(doc)
  local has = false
  doc.blocks:walk { Div = function(d) if d.identifier == 'si-contents' then has = true end end }
  if not has then return end
  -- sections: level-1 headings whose text starts with a section number ("S4." or "4.")
  local secs = {}
  local cur = nil
  local nfig, ntbl = 0, 0
  doc.blocks:walk {
    traverse = 'topdown',
    Header = function(h)
      if h.level == 1 then
        local txt = stringify(h.content)
        if txt:match('^S?%d+%.%s') then
          cur = { title = txt, fig = {}, tbl = {} }
          secs[#secs + 1] = cur
        else
          cur = nil
        end
      end
    end,
    Figure = function(f)
      nfig = nfig + 1
      if cur then cur.fig[#cur.fig + 1] = nfig end
      return nil, false
    end,
    Table = function(t)
      if stringify(t.caption.long or {}) == '' then return nil, false end
      ntbl = ntbl + 1
      if cur then cur.tbl[#cur.tbl + 1] = ntbl end
      return nil, false
    end,
  }
  local paras = pandoc.Blocks {}
  for _, sc in ipairs(secs) do
    local parts = {}
    local t = float_list('tbl', sc.tbl)
    local f = float_list('fig', sc.fig)
    if t then parts[#parts + 1] = t end
    if f then parts[#parts + 1] = f end
    local inl = pandoc.Inlines(text_to_inlines(sc.title))
    if #parts > 0 then
      inl[#inl] = pandoc.Str(stringify(inl[#inl]) .. ':')
      inl:insert(pandoc.Space())
      inl:extend(text_to_inlines(table.concat(parts, '; ')))
    end
    paras:insert(pandoc.Para(inl))
  end
  doc.blocks = doc.blocks:walk {
    Div = function(d)
      if d.identifier == 'si-contents' then
        return pandoc.Div(paras, pandoc.Attr('si-contents', { 'contents' }, { ['custom-style'] = 'Contents List' }))
      end
    end,
  }
end

---------------------------------------------------------------------------------------------------
local function json_escape(s)
  return '"' .. s:gsub('[%c"\\]', function(c) return ('\\u%04x'):format(c:byte()) end) .. '"'
end

function Pandoc(doc)
  if doc.meta['number-prefix'] then PREFIX = stringify(doc.meta['number-prefix']) end
  if doc.meta['xref-style'] and stringify(doc.meta['xref-style']) == 'acs' then
    NAMES.fig = { 'Figure', 'Figures' }
    NAMES.eq = { 'eq', 'eqs' }
    FIG_LABEL = 'Figure'
    EQ_PAREN = false
  end
  if doc.meta['max-figure-width'] then max_w_cm = tonumber(stringify(doc.meta['max-figure-width'])) or max_w_cm end
  if doc.meta['max-figure-height'] then max_h_cm = tonumber(stringify(doc.meta['max-figure-height'])) or max_h_cm end
  if IS_LATEX and doc.meta['pdf-max-figure-height'] then
    max_h_cm = math.min(max_h_cm, tonumber(stringify(doc.meta['pdf-max-figure-height'])) or max_h_cm)
  end

  collect(doc.blocks)
  contents_list(doc)
  doc = front_matter(doc)
  if IS_LATEX and PREFIX ~= '' then
    doc.blocks:insert(1, pandoc.RawBlock('latex', ('\\renewcommand{\\thefigure}{%s\\arabic{figure}}'
      .. '\\renewcommand{\\thetable}{%s\\arabic{table}}'):format(PREFIX, PREFIX)))
  end
  doc.blocks = doc.blocks:walk {
    traverse = 'topdown',
    Figure = function(f) return Figure(f), false end,   -- caption contents are walked below
  }
  -- the other transformations (typewise walk; Figure handled above so each figure is counted once)
  doc.blocks = doc.blocks:walk { Table = Table }
  if IS_LATEX then
    doc.blocks = join_display_paras(doc.blocks) or doc.blocks
    doc.blocks = doc.blocks:walk { Blocks = join_display_paras }
  end
  doc.blocks = doc.blocks:walk { Para = equation_para, Plain = function(p)
    if #find_equations(p.content) > 0 then err('display equation inside a table cell or list item is not supported') end
  end }
  doc.blocks = doc.blocks:walk { Cite = Cite }
  doc.blocks = doc.blocks:walk { Str = Str, Header = Header }
  doc.blocks = doc.blocks:walk { Image = function(img) return size_image(img) end }
  if doc.meta.title and pandoc.utils.type(doc.meta.title) == 'Inlines' then
    doc.meta.title = doc.meta.title:walk { Cite = Cite, Str = Str }
  end

  for _, w in ipairs(warnings) do io.stderr:write('[crossref] warning: ' .. w .. '\n') end
  local mapfile = doc.meta['crossref-map'] and stringify(doc.meta['crossref-map'])
  if mapfile and mapfile ~= '' then
    local f = io.open(mapfile, 'w')
    if f then
      local rows = {}
      for _, l in ipairs(order) do
        rows[#rows + 1] = ('  %s: {"kind": %s, "num": %s}'):format(json_escape(l), json_escape(labels[l].kind),
          json_escape(labels[l].num))
      end
      f:write('{\n' .. table.concat(rows, ',\n') .. '\n}\n')
      f:close()
    end
  end
  if #errors > 0 then
    for _, e in ipairs(errors) do io.stderr:write('[crossref] ERROR: ' .. e .. '\n') end
    error(('crossref.lua: %d error(s); see messages above'):format(#errors))
  end
  return doc
end
