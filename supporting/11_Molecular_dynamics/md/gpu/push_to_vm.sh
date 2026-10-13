#!/bin/bash
# push_to_vm.sh (run on the Pi): copy to the rented GPU machine exactly what it needs, to the SAME ABSOLUTE PATHS,
# with a checksum list that the machine verifies. Part of md/gpu/GPU_RUN.md, section 2.
#
#   md/gpu/push_to_vm.sh                      DRY RUN (the default): what would be copied and how much, measured
#   md/gpu/push_to_vm.sh --dry-run            locally with rsync's dry run. No machine is contacted, nothing is written.
#   md/gpu/push_to_vm.sh --preview [--ip A.B.C.D]   ASK THE MACHINE what --yes would do there, and change nothing:
#                                             for every part, which files would be NEW there, which would REPLACE a
#                                             file of other content, which would only get the Pi's time stamp (the
#                                             content is the same), and which stay because they are newer on the
#                                             machine. rsync dry runs with checksums (-n -c, with and without
#                                             --update): nothing is copied, written, installed or logged on either side.
#                                             With --list every name is printed. LOOK AT THIS BEFORE A PUSH TO A
#                                             MACHINE THAT IS WORKING: --yes replaces scripts and inputs under running tasks.
#   md/gpu/push_to_vm.sh --yes [--ip A.B.C.D] copy, send the checksum list, verify it on the machine
#   md/gpu/push_to_vm.sh --verify [--ip ...]  only verify the last checksum list on the machine again
#   md/gpu/push_to_vm.sh --rehearse           make the checksum list and check it against the Pi's own files; no machine
#                                             is contacted and nothing is copied (a test of the list)
#
# Options: --ip ADDRESS          address of the machine; default: key "ip" of /root/.config/cqc-sfa/vm_gpu.json, the state
#                                file of ops/ibm/cqc_gpu.py (its key "ssh_user" names the login, default ubuntu)
#          --include-pi-state    also copy the Pi's bookkeeping and small-model files of part D (see below). Only for
#                                moving such a job to the machine for good; it OVERWRITES the machine's own files.
#          --list                with the dry run: print every file
#
# WHAT IS COPIED (nothing is ever deleted on the machine; a later call copies only what changed ON THE PI: a file
# that is NEWER on the machine than on the Pi is left as it is there (rsync --update), named in
# push_kept_newer.txt and left out of the checksum list. The machine's own tasks rewrite files inside md/, for
# example md/qm/results/CA2_qm_summary.json; a later push must not put the Pi's older ones over them)
#  A  <repository>/md/              scripts, protocol, ligand and system inputs, the files of MCPB.py step 1, md/gpu/
#                                   (without __pycache__)
#  B  <repository>/ops/vm/setup_gpu.sh
#  C  /root/cqc-sfa-md/systems/     the built systems. NOT: build/ (work files of tleap), _checks/ (test runs),
#                                   _selftest_bonded_* (mocks marked MOCK_NOT_A_MODEL, which must never be simulated)
#  D  /root/cqc-sfa-md/qm/          checkpoints and results of the quantum chemistry. NOT: PySCF scratch (tmp* and
#                                   every file without an extension: the 17 GB of fitting integrals, the 2 GB
#                                   file of a gradient), and NOT what the machine's own jobs write
#                                   under the same names: the bookkeeping (status.json, *.pid, run_qm.log, *.stdout,
#                                   queue.log, STOP) and the files of the small-model jobs (prod/opt*, prod/hess*).
#                                   A later call must not overwrite the machine's live files with the Pi's
#  E  /root/cqc-sfa-md/ligand/, /root/cqc-sfa-md/system/    small work files that md/ligand_common.py and
#                                   md/check_tleap.py expect at these paths
#
# CHECKSUMS. Before copying, sha256 of every file goes to /root/cqc-sfa-md/status/push_manifest.sha256 (absolute
# paths, the format of sha256sum). Files that changed on the Pi while the copy ran (logs and checkpoints of running
# jobs) are taken out of the list and named in push_volatile.txt, and so are the files that stayed as they are on
# the machine because they are newer there (push_kept_newer.txt). The machine then runs `sha256sum -c` on the list.
#
# SSH: the project's key /root/.ssh/aitf_pi and options, user ubuntu, rsync on the far side through `sudo -n` (the
# target folders belong to root; no sudo if the login is root). The script prints no secrets, contacts nothing but the machine, and stops at the
# first error. Log (not for the dry run): /root/cqc-sfa-md/status/push_to_vm.log, times in UTC and IST.
# PUSH_STATUS_DIR names another folder on the Pi for the checksum list and the log (tests).
set -euo pipefail

REPO=/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate
BULK=/root/cqc-sfa-md
STATE=/root/.config/cqc-sfa/vm_gpu.json
KEY=/root/.ssh/aitf_pi
VM_USER=""                                       # from the state file; ubuntu if it names none
SSH_OPTS=(-i "$KEY" -o IdentitiesOnly=yes -o BatchMode=yes -o ConnectTimeout=15 -o StrictHostKeyChecking=accept-new -o ServerAliveInterval=30)
STATUS=${PUSH_STATUS_DIR:-$BULK/status}          # where the checksum list and the log go (on the Pi)
MODE=dry; IP=""; PI_STATE=0; LIST=0

while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) MODE=dry ;;
    --yes) MODE=push ;;
    --preview) MODE=preview ;;
    --verify) MODE=verify ;;
    --rehearse) MODE=rehearse ;;
    --ip) IP=${2:?--ip needs an address}; shift ;;
    --include-pi-state) PI_STATE=1 ;;
    --list) LIST=1 ;;
    -h|--help) sed -n '2,/^set -euo pipefail/p' "$0" | sed '$d'; exit 0 ;;
    *) echo "unknown option $1 (see --help)" >&2; exit 2 ;;
  esac
  shift
done

log() { echo "[$(date -u '+%Y-%m-%dT%H:%M:%SZ') | $(date -u -d '+5 hours 30 minutes' '+%-d %b %Y, %-I:%M:%S %p IST')] $*"; }

# the parts: source | destination | rsync filters
PARTS=("A md"       "B setup_gpu.sh" "C systems"  "D qm"  "E ligand"  "E system")
SRC=("$REPO/md/"    "$REPO/ops/vm/setup_gpu.sh" "$BULK/systems/" "$BULK/qm/" "$BULK/ligand/" "$BULK/system/")
DST=("$REPO/md/"    "$REPO/ops/vm/setup_gpu.sh" "$BULK/systems/" "$BULK/qm/" "$BULK/ligand/" "$BULK/system/")
filters() {   # filters <index>: sets the array FILT, the rsync arguments that select the files of a part
  case "$1" in
    0) FILT=('--exclude=__pycache__/' '--exclude=*.pyc') ;;
    2) FILT=('--exclude=/_checks/' '--exclude=/_selftest_bonded_mock/' '--exclude=/_selftest_bonded_synthetic/' '--exclude=build/') ;;
    3) FILT=('--exclude=tmp*')
       if [ "$PI_STATE" != 1 ]; then
         FILT+=('--exclude=status.json' '--exclude=*.pid' '--exclude=run_qm.log' '--exclude=*.stdout' '--exclude=queue.log' '--exclude=STOP'
                '--exclude=/prod/opt*' '--exclude=/prod/hess*')
       fi
       # PySCF's scratch files have no extension and not always the prefix tmp (seen: t71ug9bl, 2 GB, during a
       # gradient): of this folder only files with an extension are copied
       FILT+=('--include=*/' '--include=*.*' '--exclude=*') ;;
    *) FILT=() ;;
  esac
}

EMPTY=$(mktemp -d)                      # an empty folder: the local dry run "copies" into it
TMP=$(mktemp -d)
trap 'rm -rf "$EMPTY" "$TMP"' EXIT

file_list() {   # file_list <index>: absolute paths of the regular files of a part, one per line
  local i=$1 src=${SRC[$1]}
  if [ -f "$src" ]; then echo "$src"; return; fi
  filters "$i"
  rsync -a -n --out-format='%n' "${FILT[@]}" "$src" "$EMPTY/" | while IFS= read -r name; do
    [ -f "$src$name" ] && [ ! -L "$src$name" ] && echo "$src$name"
  done
  return 0
}

measure() {   # the dry run: files and bytes of every part, measured locally
  local i n bytes total_n=0 total_b=0
  printf '%-16s %8s %12s   %s\n' part files size source
  for i in "${!PARTS[@]}"; do
    [ -e "${SRC[$i]}" ] || { echo "MISSING: ${SRC[$i]}" >&2; exit 1; }
    file_list "$i" > "$TMP/list_$i"
    n=$(wc -l < "$TMP/list_$i")
    bytes=$(tr '\n' '\0' < "$TMP/list_$i" | du -cb --files0-from=- | tail -n 1 | cut -f1)
    total_n=$((total_n + n)); total_b=$((total_b + bytes))
    printf '%-16s %8d %9.1f MB   %s\n' "${PARTS[$i]}" "$n" "$(echo "$bytes" | awk '{print $1 / 1e6}')" "${SRC[$i]}"
    if [ "$LIST" = 1 ]; then sed 's/^/      /' "$TMP/list_$i"; fi
  done
  printf '%-16s %8d %9.1f MB\n' TOTAL "$total_n" "$(echo "$total_b" | awk '{print $1 / 1e6}')"
  echo "largest files:"
  cat "$TMP"/list_* | tr '\n' '\0' | du -b --files0-from=- | sort -rn | awk 'NR <= 5 {printf "  %8.1f MB  %s\n", $1 / 1e6, $2}'
  echo "left out on purpose (sizes on the Pi): $(du -sh "$BULK/qm" 2> /dev/null | cut -f1) in $BULK/qm before the filters," \
       "of which PySCF scratch (tmp*, files without an extension): $(find "$BULK/qm" -type f \( -name 'tmp*' -o ! -name '*.*' \) -printf '%s\n' 2> /dev/null | awk '{s += $1} END {printf "%.1f GB", s / 1e9}')"
}

address() {
  if [ -s "$STATE" ]; then
    [ -n "$IP" ] || IP=$(python3 -c "import json, sys; print(json.load(open(sys.argv[1])).get('ip') or '')" "$STATE")
    VM_USER=$(python3 -c "import json, sys; print(json.load(open(sys.argv[1])).get('ssh_user') or '')" "$STATE")
  fi
  VM_USER=${VM_USER:-ubuntu}
  [ -n "$IP" ] || { echo "no address: give --ip, or let ops/ibm/cqc_gpu.py write \"ip\" into $STATE" >&2; exit 2; }
  case "$IP" in *[!0-9.]*) echo "not an IPv4 address: $IP" >&2; exit 2 ;; esac
  case "$VM_USER" in *[!a-z0-9_-]*) echo "odd login name in $STATE" >&2; exit 2 ;; esac
  if [ "$VM_USER" = root ]; then SUDO=""; RSYNC_PATH="rsync"; else SUDO="sudo -n "; RSYNC_PATH="sudo -n rsync"; fi
}
remote() { ssh "${SSH_OPTS[@]}" "$VM_USER@$IP" "$@"; }

verify_remote() {
  log "verifying on the machine: sha256sum -c /root/cqc-sfa-md/status/push_manifest.sha256"
  if remote "${SUDO}sh -c 'sha256sum --quiet -c $BULK/status/push_manifest.sha256 && wc -l < $BULK/status/push_manifest.sha256'"; then
    log "VERIFIED: every file of the list has its checksum on the machine"
  else
    log "VERIFICATION FAILED: see the lines above; copy again (--yes) and look at $STATUS/push_volatile.txt"; exit 1
  fi
}

make_manifest() {   # checksums of every file of the lists, and size and time of each as they are now
  log "checksums of the files on the Pi"
  : > "$TMP/manifest"; : > "$TMP/stat_before"
  local i
  for i in "${!PARTS[@]}"; do
    tr '\n' '\0' < "$TMP/list_$i" | xargs -0 -r nice -n 10 sha256sum >> "$TMP/manifest"
    tr '\n' '\0' < "$TMP/list_$i" | xargs -0 -r stat -c '%n|%s|%y' >> "$TMP/stat_before"
  done
}

finish_manifest() {   # files that changed on the Pi since their checksum was taken cannot be verified: out of the list, and named
  local i
  : > "$TMP/stat_after"
  for i in "${!PARTS[@]}"; do
    tr '\n' '\0' < "$TMP/list_$i" | xargs -0 -r stat -c '%n|%s|%y' >> "$TMP/stat_after" 2> /dev/null || true
  done
  python3 - "$TMP/manifest" "$TMP/stat_before" "$TMP/stat_after" "$STATUS/push_manifest.sha256" "$STATUS/push_volatile.txt" \
            "$TMP/kept_newer" "$STATUS/push_kept_newer.txt" <<'EOF'
import os, sys
manifest, before, after, out, volatile, kept_in, kept_out = sys.argv[1:8]
b = dict(l.rstrip("\n").split("|", 1) for l in open(before))
a = dict(l.rstrip("\n").split("|", 1) for l in open(after))
changed = sorted(p for p in b if a.get(p) != b[p])
kept = sorted({l.rstrip("\n") for l in open(kept_in) if l.strip()}) if os.path.exists(kept_in) else []
drop = set(changed) | set(kept)
keep = [l for l in open(manifest) if l.rstrip("\n").split("  ", 1)[1] not in drop]
open(out + ".tmp", "w").writelines(keep)
os.replace(out + ".tmp", out)
open(volatile, "w").write("".join(p + "\n" for p in changed))
open(kept_out, "w").write("".join(p + "\n" for p in kept))
print(f"    checksum list: {len(keep)} files; changed on the Pi during the copy and left out: {len(changed)}"
      + (" (" + ", ".join(os.path.basename(p) for p in changed[:6]) + ")" if changed else ""))
if kept:
    print(f"    NEWER ON THE MACHINE, left as they are there and out of the list: {len(kept)} (named in {kept_out})")
    for p in kept[:12]:
        print(f"      {p}")
EOF
}

if [ "$MODE" = dry ]; then
  echo "DRY RUN (nothing is copied, no machine is contacted). With --yes these files go to the same paths on the machine:"
  measure
  echo "to copy: $0 --yes [--ip ADDRESS]"
  exit 0
fi

if [ "$MODE" != rehearse ]; then address; fi      # a call without an address ends here and writes nothing

if [ "$MODE" = preview ]; then
  # The machine is asked, nothing is changed: rsync -n on both sides, no mkdir, no install, no checksum list, no log.
  log "PREVIEW of a push to $VM_USER@$IP: the machine is asked; nothing is copied or written there."
  remote "command -v rsync > /dev/null" ||
    { echo "no preview: the machine does not answer on SSH (address $IP, key $KEY, user $VM_USER), or rsync is not installed there (--yes installs it)" >&2; exit 1; }
  t_new=0; t_rep=0; t_time=0; t_kept=0; pv_rc=0
  for i in "${!PARTS[@]}"; do
    [ -e "${SRC[$i]}" ] || { echo "MISSING: ${SRC[$i]}" >&2; exit 1; }
    filters "$i"
    ok=1
    for u in "" --update; do      # as --yes does it: what would go without --update and not with it is newer on the machine
      rsync -a -n -c $u --out-format='%i|%n' -e "ssh ${SSH_OPTS[*]}" --rsync-path="$RSYNC_PATH" "${FILT[@]}" "${SRC[$i]}" "$VM_USER@$IP:${DST[$i]}" \
        > "$TMP/pv${u:+_update}" 2> "$TMP/pv_err" || ok=0
    done
    if [ "$ok" = 0 ]; then
      echo "part ${PARTS[$i]}: the machine could not be asked: $(tail -n 1 "$TMP/pv_err")"; pv_rc=1; continue
    fi
    awk -F'|' '$1 ~ /^<f\+/ {print $2}' "$TMP/pv_update" | sort > "$TMP/pv_new"
    awk -F'|' '$1 ~ /^<f/ && $1 !~ /^<f\+/ {print $2}' "$TMP/pv_update" | sort > "$TMP/pv_rep"
    awk -F'|' '$1 ~ /^\.f/ && substr($1, 5, 1) == "t" {print $2}' "$TMP/pv_update" | sort > "$TMP/pv_time"
    awk -F'|' '$1 ~ /^<f/ && $1 !~ /^<f\+/ {print $2}' "$TMP/pv" | sort | comm -23 - "$TMP/pv_rep" > "$TMP/pv_kept"
    n_new=$(wc -l < "$TMP/pv_new"); n_rep=$(wc -l < "$TMP/pv_rep"); n_time=$(wc -l < "$TMP/pv_time"); n_kept=$(wc -l < "$TMP/pv_kept")
    t_new=$((t_new + n_new)); t_rep=$((t_rep + n_rep)); t_time=$((t_time + n_time)); t_kept=$((t_kept + n_kept))
    echo "part ${PARTS[$i]} (${DST[$i]}): $n_new new on the machine, $n_rep would REPLACE a file of other content, $n_time time stamp only, $n_kept newer on the machine (stay)"
    sed 's/^/    REPLACE  /' "$TMP/pv_rep"
    sed 's/^/    stays (newer on the machine, other content)  /' "$TMP/pv_kept"
    if [ "$LIST" = 1 ]; then
      sed 's/^/    new  /' "$TMP/pv_new"; sed 's/^/    time stamp only  /' "$TMP/pv_time"
    else
      head -n 8 "$TMP/pv_new" | sed 's/^/    new  /'
      [ "$n_new" -le 8 ] || echo "    new  ... and $((n_new - 8)) more (--list prints every name)"
    fi
  done
  echo "TOTAL: $t_new new, $t_rep would REPLACE a file of other content on the machine, $t_time time stamp only, $t_kept stay as they are there."
  echo "Nothing was changed. With --yes the files marked REPLACE are put over the machine's while its tasks run: a script or an"
  echo "input that a task reads at its next step is then another one. Choose the moment (md/gpu/GPU_RUN.md, section 2)."
  exit "$pv_rc"
fi
mkdir -p "$STATUS"
exec > >(tee -a "$STATUS/push_to_vm.log") 2>&1

if [ "$MODE" = rehearse ]; then
  log "push_to_vm.sh rehearsal: the checksum list is made and checked against the Pi's own files; no machine is contacted, nothing is copied"
  measure
  make_manifest
  finish_manifest
  if (cd / && sha256sum --quiet -c "$STATUS/push_manifest.sha256"); then
    log "REHEARSAL PASSED: $(wc -l < "$STATUS/push_manifest.sha256") files of $STATUS/push_manifest.sha256 have their checksums (changed meanwhile and left out: $(wc -l < "$STATUS/push_volatile.txt"))"
    exit 0
  fi
  log "REHEARSAL FAILED: see the lines above"; exit 1
fi

log "push_to_vm.sh $MODE to $VM_USER@$IP"
remote true || { log "the machine does not answer on SSH (address $IP, key $KEY, user $VM_USER); nothing was copied"; exit 1; }
if [ "$MODE" = verify ]; then verify_remote; exit 0; fi

measure
make_manifest

log "making sure rsync and the target folders exist on the machine"
remote "command -v rsync > /dev/null || (${SUDO}apt-get -o DPkg::Lock::Timeout=300 update -qq && ${SUDO}env DEBIAN_FRONTEND=noninteractive apt-get -o DPkg::Lock::Timeout=300 install -y -qq rsync)"
remote "${SUDO}mkdir -p $REPO/md $REPO/ops/vm $BULK/systems $BULK/qm $BULK/ligand $BULK/system $BULK/status $BULK/runs $BULK/logs"

: > "$TMP/kept_newer"
for i in "${!PARTS[@]}"; do
  log "copying part ${PARTS[$i]}: ${SRC[$i]} -> $IP:${DST[$i]}"
  filters "$i"
  # A file that is newer on the machine than on the Pi stays as it is there (--update below). Which files those are
  # is found first, by two dry runs that differ only in --update: what the first would send and the second would not.
  for u in "" --update; do
    rsync -a -n $u --out-format='%n' -e "ssh ${SSH_OPTS[*]}" --rsync-path="$RSYNC_PATH" "${FILT[@]}" "${SRC[$i]}" "$VM_USER@$IP:${DST[$i]}" |
      { grep -v '/$' || true; } | sort > "$TMP/would${u:+_update}"
  done
  if [ -f "${SRC[$i]}" ]; then pre=$(dirname "${SRC[$i]}")/; else pre=${SRC[$i]}; fi
  comm -23 "$TMP/would" "$TMP/would_update" | sed "s|^|$pre|" >> "$TMP/kept_newer"
  rsync -a --update --partial --stats -e "ssh ${SSH_OPTS[*]}" --rsync-path="$RSYNC_PATH" "${FILT[@]}" "${SRC[$i]}" "$VM_USER@$IP:${DST[$i]}" |
    { grep -E 'Number of regular files transferred|Total transferred file size|Total bytes sent' || true; } | sed 's/^/    /'
done

finish_manifest
rsync -a -e "ssh ${SSH_OPTS[*]}" --rsync-path="$RSYNC_PATH" "$STATUS/push_manifest.sha256" "$STATUS/push_volatile.txt" "$STATUS/push_kept_newer.txt" "$VM_USER@$IP:$BULK/status/"
verify_remote
log "PUSH COMPLETE to $IP. Next on the machine: ${SUDO}bash $REPO/ops/vm/setup_gpu.sh (md/gpu/GPU_RUN.md, section 2)"
