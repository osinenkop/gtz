# D6/S6 additional-prime support batch

Started: 2026-08-10 08:38 UTC.

Purpose: add CRT modulus for the high-value support families whose
three-prime rational reconstruction failed:

- size-7 `78612`, support candidate `19`, extra primes `32029,32051,32057`
- size-8 `79656`, support candidate `6`, extra primes `32029,32051,32057`

Watched jobs are listed in `jobs.txt`.  The watcher writes:

- `latest_summary.md`
- `state.json`
- `raw/<job>/job_get_latest.txt`
- `logs/<job>.log`

Successful artifacts are downloaded to `code/mlcore/downloads/<job>/GTZ_OUT.zip`
and extracted into `code/sage/out`.
