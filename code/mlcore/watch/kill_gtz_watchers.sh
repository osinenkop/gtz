#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/p.osinenko/Documents/gtz}"
DRY_RUN=0
PATTERN='code/mlcore/watch/(mlcore_watch|d6_followup_loop)[.]py'
DIR_FILTER=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --followup-only)
      PATTERN='code/mlcore/watch/d6_followup_loop[.]py'
      shift
      ;;
    --watchers-only)
      PATTERN='code/mlcore/watch/mlcore_watch[.]py'
      shift
      ;;
    --dir)
      if [[ "$#" -lt 2 ]]; then
        echo "--dir requires a directory name or path" >&2
        exit 2
      fi
      if [[ "$2" = /* ]]; then
        DIR_FILTER="$2"
      else
        DIR_FILTER="$ROOT/code/mlcore/watch/$2"
      fi
      shift 2
      ;;
    *)
      echo "unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

mapfile -t candidates < <(
  ps -eo pid=,pgid=,args= |
    awk -v root="$ROOT" -v pattern="$PATTERN" -v dir_filter="$DIR_FILTER" '
      $0 ~ pattern && index($0, "--root " root) {
        if (dir_filter != "" && index($0, "--watch-dir " dir_filter) == 0 && index($0, "--run-dir " dir_filter) == 0) {
          next
        }
        pid=$1; pgid=$2;
        $1=""; $2="";
        sub(/^  */, "", $0);
        print pid " " pgid " " $0;
      }
    '
)

if [[ "${#candidates[@]}" -eq 0 ]]; then
  echo "No GTZ MLCore watcher processes found for root: $ROOT"
  exit 0
fi

printf '%s\n' "${candidates[@]}"

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "Dry run only; no processes killed."
  exit 0
fi

for row in "${candidates[@]}"; do
  pid="${row%% *}"
  rest="${row#* }"
  pgid="${rest%% *}"

  if [[ "$pgid" =~ ^[0-9]+$ ]] && [[ "$pgid" != "0" ]]; then
    kill -TERM "-$pgid" 2>/dev/null || kill -TERM "$pid" 2>/dev/null || true
  else
    kill -TERM "$pid" 2>/dev/null || true
  fi
done

sleep 2

still_running=0
for row in "${candidates[@]}"; do
  pid="${row%% *}"
  if ps -p "$pid" >/dev/null 2>&1; then
    still_running=1
    echo "Still running after TERM: $pid"
  fi
done

if [[ "$still_running" -eq 0 ]]; then
  echo "Stopped GTZ MLCore watcher processes."
fi
