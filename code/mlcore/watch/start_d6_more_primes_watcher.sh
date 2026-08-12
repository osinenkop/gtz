#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/p.osinenko/Documents/gtz}"
RUN_DIR="$ROOT/code/mlcore/watch/d6_more_primes"
JOBS_FILE="$RUN_DIR/jobs.txt"
PID_FILE="$RUN_DIR/watcher.pid"
STDOUT_LOG="$RUN_DIR/watcher_stdout.log"
START_LOG="$RUN_DIR/watcher_start.log"

cd "$ROOT"
mkdir -p "$RUN_DIR"

if [[ -s "$PID_FILE" ]]; then
  old_pid="$(cat "$PID_FILE")"
  if ps -p "$old_pid" >/dev/null 2>&1; then
    echo "watcher already running: $old_pid" | tee "$START_LOG"
    exit 0
  fi
fi

mapfile -t jobs < "$JOBS_FILE"
if [[ "${#jobs[@]}" -eq 0 ]]; then
  echo "no jobs listed in $JOBS_FILE" | tee "$START_LOG"
  exit 1
fi

nohup setsid python3 code/mlcore/watch/mlcore_watch.py \
  --project aida \
  --jobs "${jobs[@]}" \
  --interval-seconds 600 \
  --hours 8 \
  --root "$ROOT" \
  --watch-dir "$RUN_DIR" \
  --continue-after-breakthrough \
  --continue-until-deadline \
  --command-timeout 180 \
  --download-timeout 600 \
  > "$STDOUT_LOG" 2>&1 < /dev/null &

pid="$!"
echo "$pid" > "$PID_FILE"
sleep 2

if ps -p "$pid" >/dev/null 2>&1; then
  {
    echo "watcher started: $pid"
    ps -o pid,ppid,stat,etime,cmd -p "$pid"
  } | tee "$START_LOG"
else
  {
    echo "watcher failed to stay running: $pid"
    tail -80 "$STDOUT_LOG" || true
  } | tee "$START_LOG"
  exit 1
fi
