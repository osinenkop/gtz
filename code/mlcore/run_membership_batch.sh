#!/usr/bin/env bash
set -euo pipefail

# Run from an MLCore job image that provides Sage as `sage`.
# The expected input layout is:
#   /work/gtz/sage/{*.py,out/refine_*.json}

cd "${GTZ_SAGE_DIR:-/work/gtz/sage}"
mkdir -p out

DEGREE="${DEGREE:-5}"
TARGET_POWER="${TARGET_POWER:-1}"
SEPARATOR="${SEPARATOR:-chart-d}"
SEPARATOR_POWER="${SEPARATOR_POWER:-1}"
CHARACTERISTIC="${CHARACTERISTIC:-32003}"
CASES="${CASES:-s7,s8}"

run_case() {
  local case_id="$1"
  local input rows target out

  case "$case_id" in
    s7)
      input="out/refine_overtie_s7_78612_t2_7_10_19_p1400.json"
      rows="0,2,3,4,5,6,7,8,9,10"
      target="1"
      out="out/membership_s7_78612_omit_active1_${SEPARATOR//-/_}${SEPARATOR_POWER}_pow${TARGET_POWER}_deg${DEGREE}_p${CHARACTERISTIC}.json"
      ;;
    s8)
      input="out/refine_overtie_s8_79656_t10_11_15_p1400.json"
      rows="1,2,3,4,5,6,7,8,9,10"
      target="0"
      out="out/membership_s8_79656_omit_active0_${SEPARATOR//-/_}${SEPARATOR_POWER}_pow${TARGET_POWER}_deg${DEGREE}_p${CHARACTERISTIC}.json"
      ;;
    *)
      echo "unknown case: $case_id" >&2
      return 2
      ;;
  esac

  echo "running case=$case_id degree=$DEGREE target_power=$TARGET_POWER separator=$SEPARATOR^$SEPARATOR_POWER characteristic=$CHARACTERISTIC"
  sage -python probe_bounded_ideal_membership.py \
    --input "$input" \
    --rows "$rows" \
    --target-row "$target" \
    --target-power "$TARGET_POWER" \
    --separator "$SEPARATOR" \
    --separator-power "$SEPARATOR_POWER" \
    --degree "$DEGREE" \
    --characteristic "$CHARACTERISTIC" \
    --out "$out"
}

IFS=',' read -r -a case_list <<< "$CASES"
for case_id in "${case_list[@]}"; do
  run_case "$case_id"
done
