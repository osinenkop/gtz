# D6/S6 support-restricted separator batch

Started: 2026-08-09 22:00 UTC.

Purpose: test whether selected nonzero `F_32003` D6/S6 separator supports
persist over `F_32009` and `F_32027`, so that a CRT/rational reconstruction can
be attempted for an actual characteristic-zero local denominator.

Watched jobs are listed in `jobs.txt`.  Current status is written by
`code/mlcore/watch/mlcore_watch.py` to:

- `latest_summary.md`
- `state.json`
- `raw/<job>/job_get_latest.txt`
- `logs/<job>.log`

Selected supports:

- size-7 `78612`: candidate indices `19` and `8` from
  `code/sage/out/local_separator_s7_78612_omit_active1_D6_S6_32003.json`
- size-8 `79656`: candidate indices `6` and `8` from
  `code/sage/out/local_separator_s8_79656_omit_active0_D6_S6_32003.json`

The watcher runs for 8 hours, polls every 10 minutes, downloads successful
artifacts to `code/mlcore/downloads/<job>/GTZ_OUT.zip`, and extracts
`local_separator*.json` outputs into `code/sage/out`.
