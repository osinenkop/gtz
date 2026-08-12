#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/p.osinenko/Documents/gtz
WATCH_DIR="$ROOT/code/mlcore/watch"
PID_FILE="$WATCH_DIR/watcher.pid"
STDOUT_LOG="$WATCH_DIR/watcher_stdout.log"
START_LOG="$WATCH_DIR/watcher_start.log"

cd "$ROOT"
mkdir -p "$WATCH_DIR"

if [[ -s "$PID_FILE" ]]; then
  old_pid="$(cat "$PID_FILE")"
  if ps -p "$old_pid" >/dev/null 2>&1; then
    echo "watcher already running: $old_pid" | tee "$START_LOG"
    exit 0
  fi
fi

nohup python3 code/mlcore/watch/mlcore_watch.py \
  --project aida \
  --interval-seconds 600 \
  --hours 12 \
  --root "$ROOT" \
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
