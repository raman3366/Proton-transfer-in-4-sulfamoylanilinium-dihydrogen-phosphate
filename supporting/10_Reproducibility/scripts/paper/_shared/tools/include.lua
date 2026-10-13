--[[ include.lua: replace a line  <!-- include: tables/tab_geometry.md -->  with the parsed contents of that file.
Paths are relative to the manuscript's directory (pandoc runs there). Runs before crossref.lua, so included
tables and figures are numbered like any other. A missing file stops the build.
]]
function RawBlock(el)
  if el.format ~= 'html' then return nil end
  local path = el.text:match('^<!%-%-%s*include:%s*(.-)%s*%-%->$')
  if not path then return nil end
  local fh = io.open(path, 'r')
  if not fh then
    io.stderr:write('[crossref] ERROR: include file not found: ' .. path .. '\n')
    error('include file not found: ' .. path)
  end
  local text = fh:read('a')
  fh:close()
  return pandoc.read(text, 'markdown').blocks
end
