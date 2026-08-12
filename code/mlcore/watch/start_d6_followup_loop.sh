#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/p.osinenko/Documents/gtz}"
RUN_DIR="$ROOT/code/mlcore/watch/d6_followup"
PID_FILE="$RUN_DIR/followup.pid"
STDOUT_LOG="$RUN_DIR/followup_stdout.log"
START_LOG="$RUN_DIR/followup_start.log"
HOURS="${HOURS:-8.5}"

cd "$ROOT"
mkdir -p "$RUN_DIR"

if [[ -s "$PID_FILE" ]]; then
  old_pid="$(cat "$PID_FILE")"
  old_cmd="$(ps -p "$old_pid" -o args= 2>/dev/null || true)"
  if [[ "$old_cmd" == *"code/mlcore/watch/d6_followup_loop.py"* ]] && [[ "$old_cmd" == *"--run-dir $RUN_DIR"* ]]; then
    echo "follow-up loop already running: $old_pid" | tee "$START_LOG"
    exit 0
  fi
fi

nohup setsid python3 code/mlcore/watch/d6_followup_loop.py \
  --root "$ROOT" \
  --run-dir "$RUN_DIR" \
  --interval-seconds 600 \
  --hours "$HOURS" \
  > "$STDOUT_LOG" 2>&1 < /dev/null &

pid="$!"
echo "$pid" > "$PID_FILE"
sleep 2

if ps -p "$pid" >/dev/null 2>&1; then
  {
    echo "follow-up loop started: $pid"
    ps -o pid,ppid,stat,etime,cmd -p "$pid"
  } | tee "$START_LOG"
else
  {
    echo "follow-up loop failed to stay running: $pid"
    tail -80 "$STDOUT_LOG" || true
  } | tee "$START_LOG"
  exit 1
fi
