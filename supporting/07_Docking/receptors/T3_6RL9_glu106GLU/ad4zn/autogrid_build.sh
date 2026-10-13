#!/usr/bin/env bash
# Build of autogrid4 for linux-aarch64 from the official source (5 Oct 2026), for dock/GLU106.md of the
# SAPH repository.  Nothing outside this prefix is touched; the compilers are the system ones.
#
#   source   https://github.com/ccsb-scripps/AutoGrid, tag v4.2.9 (commit 6d2847beaeac8ff43ca99094707fd74e3ca1ff37)
#            archive/refs/tags/v4.2.9.tar.gz, sha256 d48eb3c72eb3708016ee27c297e2b7c3ad1e3d10e3ea2b0ea0cba95d4dcfbaf7
#            (the program reports itself as "AutoGrid 4.2.7.x", main.cpp $Revision: 1.160 $, as in the
#            stored receptor_tz.glg of the reported maps)
#   method   the autotools build of the README: autoreconf -i; ../configure; make
#   two departures from a plain build, both recorded here:
#   1. default_parameters.h is generated with sh instead of csh (no csh on this machine).  The
#      generator ad4_shared/paramdat2h.csh consists of echo, egrep and sed commands only, which sh
#      runs unchanged; the header is written before make, so the csh rule of the Makefile is not run.
#   2. CXXFLAGS="-O2 -ffp-contract=off" for the binary that is installed.  On aarch64 GCC contracts
#      a*b+c into a fused multiply-add by default, which x86-64 (the machine of the reported maps)
#      does not have in its base instruction set.  AutoGrid looks its pair energies up in tables
#      indexed by (int)(r*128), so one differing last bit in r moves a pair to the next table entry.
#      With contraction off both machines do the same IEEE operations in the same order.  The build
#      with the default flags (-g -O2) is kept beside it as autogrid4_default_flags for comparison.
set -eu
P=$(cd "$(dirname "$0")" && pwd)
SRC=$P/src/AutoGrid-4.2.9
build () {  # name cxxflags
  local B=$P/build_$1
  rm -rf "$B"; mkdir -p "$B"; cd "$B"
  sh "$SRC/ad4_shared/paramdat2h.csh" "$SRC/ad4_shared/AD4_parameters.dat" "$SRC/ad4_shared/AD4.1_bound.dat" > default_parameters.h
  CXXFLAGS="$2" "$SRC/configure" --prefix="$P" > configure.log 2>&1
  touch default_parameters.h
  nice -n 10 make -j2 > make.log 2>&1
  ls -la autogrid4
}
build nofma "-O2 -ffp-contract=off"
cp "$P/build_nofma/autogrid4" "$P/bin/autogrid4"
if [ "${WITH_DEFAULT:-1}" = 1 ]; then
  build default "-g -O2"
  cp "$P/build_default/autogrid4" "$P/bin/autogrid4_default_flags"
fi
sha256sum "$P"/bin/autogrid4*
