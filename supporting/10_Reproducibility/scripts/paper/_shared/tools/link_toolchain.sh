#!/usr/bin/env bash
# link_toolchain.sh -- point this repo's git-ignored tools/ at the paper toolchain of the sister project.
#
# Why: build.py, docx_check.py, make_reference_docx.py and figkit look for their third-party tools at
# repo-relative paths (tools/bin/pandoc, tools/paper-env, tools/texlive-extracted, tools/fonts/extracted),
# exactly as in cqc-dft. That toolchain (pandoc 3.11, a Python 3.12 venv with lxml/python-docx/matplotlib/
# pyvista/rdkit, the extracted luaotfload tree and the TeX Gyre / Liberation fonts, ~1.6 GB) already exists
# on the Pi in /root/workspace/cqc-dft/tools, so it is linked, not reinstalled. tools/ is git-ignored.
#
#   bash paper/_shared/tools/link_toolchain.sh [SISTER_TOOLS_DIR]      (idempotent; default below)
#
# The sister repo is only read, never modified. On a machine without it, install the same pieces under
# tools/ instead (pandoc binary in tools/bin/, `python3 -m venv tools/paper-env` + pip install lxml
# python-docx matplotlib numpy scipy pyvista rdkit gemmi, and the texlive-luatex / fonts-texgyre /
# fonts-liberation .deb files extracted with dpkg-deb -x as in cqc-dft/tools/dl/debs).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
SRC="${1:-/root/workspace/cqc-dft/tools}"
for need in bin/pandoc paper-env/bin/python texlive-extracted fonts/extracted; do
  [[ -e "$SRC/$need" ]] || { echo "missing $SRC/$need" >&2; exit 1; }
done
mkdir -p "$REPO/tools/bin" "$REPO/tools/fonts"
link() {   # link TARGET LINKNAME: create or refresh a symlink, never overwrite a real file/dir
  if [[ -e "$2" && ! -L "$2" ]]; then echo "exists and is not a symlink, left alone: $2" >&2; return; fi
  ln -sfn "$1" "$2"; echo "$2 -> $1"
}
link "$SRC/bin/pandoc"        "$REPO/tools/bin/pandoc"
link "$SRC/paper-env"         "$REPO/tools/paper-env"
link "$SRC/texlive-extracted" "$REPO/tools/texlive-extracted"
link "$SRC/fonts/extracted"   "$REPO/tools/fonts/extracted"
"$REPO/tools/bin/pandoc" --version | head -1
"$REPO/tools/paper-env/bin/python" -c "import lxml, docx, matplotlib; print('paper-env OK')"
