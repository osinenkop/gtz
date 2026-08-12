#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/p.osinenko/Documents/gtz}"
PROJECT="${PROJECT:-aida}"
RUN_DIR="${RUN_DIR:-$ROOT/code/mlcore/watch/d6_support}"
JOBS_FILE="$RUN_DIR/jobs.txt"
SUBMIT_LOG="$RUN_DIR/submit.log"
PRESET="code/mlcore/sage_local_separator.mlcore-preset.yaml"

mkdir -p "$RUN_DIR"
touch "$JOBS_FILE" "$SUBMIT_LOG"
cd "$ROOT"

submit_support_job() {
  local gtz_case="$1"
  local prime="$2"
  local candidate="$3"
  local generator_rank="$4"
  local support_from="$5"
  local prefix="gtz-d6sup-${gtz_case}-p${prime}-c${candidate}"

  if grep -q "^${prefix}-" "$JOBS_FILE"; then
    echo "skip existing prefix ${prefix}" | tee -a "$SUBMIT_LOG"
    return
  fi

  echo "submitting ${prefix}" | tee -a "$SUBMIT_LOG"
  local output
  output="$(
    mlc job submit -d -p "$PROJECT" \
      --preset-file "$PRESET" \
      --custom-preset-context . \
      --generate-name "$prefix" \
      --env "GTZ_CASE=${gtz_case}" \
      --env "DEGREE=6" \
      --env "SEPARATOR_DEGREE=6" \
      --env "CHARACTERISTIC=${prime}" \
      --env "SUPPORT_FROM=${support_from}" \
      --env "SUPPORT_CANDIDATE_INDEX=${candidate}" \
      --env "KNOWN_GENERATOR_RANK=${generator_rank}" \
      --env "TERM_LIMIT=1200" \
      --env "CANDIDATE_LIMIT=10" \
      - 2>&1
  )"
  printf '%s\n' "$output" | tee -a "$SUBMIT_LOG"

  local job_name
  job_name="$(printf '%s\n' "$output" | awk '/Job .* created at/ {print $4}' | tail -1)"
  if [[ -z "$job_name" ]]; then
    echo "failed to parse job name for ${prefix}" | tee -a "$SUBMIT_LOG"
    return 1
  fi
  printf '%s\n' "$job_name" >> "$JOBS_FILE"
}

S7_SUPPORT="out/local_separator_s7_78612_omit_active1_D6_S6_32003.json"
S8_SUPPORT="out/local_separator_s8_79656_omit_active0_D6_S6_32003.json"

# Candidate choices:
# - c19/c6 are the largest saved |s(root)| candidates for size 7/8 respectively.
# - c8 is a small-support candidate in both cases, useful for rational reconstruction.
for prime in 32003 32009 32027; do
  submit_support_job s7 "$prime" 19 79399 "$S7_SUPPORT"
  submit_support_job s7 "$prime" 8 79399 "$S7_SUPPORT"
  submit_support_job s8 "$prime" 6 79071 "$S8_SUPPORT"
  submit_support_job s8 "$prime" 8 79071 "$S8_SUPPORT"
done

sort -u "$JOBS_FILE" -o "$JOBS_FILE"
echo "wrote $JOBS_FILE" | tee -a "$SUBMIT_LOG"
