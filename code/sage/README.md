# Sage/Singular semialgebraic workspace

This folder contains the first exact-algebra layer for the GTZ(6,3)
semialgebraic route.

Run scripts through the Sage conda environment:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/export_semialgebraic_system.py --help
```

## What is generated

`export_semialgebraic_system.py` writes:

- a JSON summary with active triples, variable counts, equation counts, and
  degrees;
- a `.sing` file that can be run by Singular.

The supported modes are:

- `det`: saturated active determinant equations
  `det(6*N_TT-d*I)/d^2=0` in the standard chart `Y=[I;Z]`, where
  `d=det(Y^T Y)`.  The saturation is valid on the full-rank chart `d!=0` and
  reduces the determinant degree from `18` to `6`; use `--expand-det` when
  asking Sage for dimensions or Groebner bases.  Use `--invert-d` to add
  `u0*d-1=0`, which removes algebraic components on the chart boundary `d=0`.
- `nonsharp`: the simple-active cofactor-patch tangent-witness system, with
  inequalities `q_T<=0` encoded as `q_T+s_T^2=0`.
- `kernel`: a lower-degree simple-active tangent-witness system with explicit
  kernel-vector variables `w`.  This has more variables than the cofactor model,
  but avoids the cofactor degree blow-up and patch branching.

The `nonsharp` and `kernel` modes are algebraic encodings/relaxations of the
simple-active obstruction.  PSD, simplicity, and inactive inequalities are not
solved here yet; this is the export/probe layer for exact backends.

## Examples

Small determinant smoke test:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/export_semialgebraic_system.py \
  --mode det --active-indices 0 --dimension \
  --out-prefix code/sage/out/smoke_det_A0
```

One-triple cofactor-patch nonsharp smoke test:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/export_semialgebraic_system.py \
  --mode nonsharp --active-indices 0 --row-pairs 01 --patch-inverses \
  --out-prefix code/sage/out/smoke_nonsharp_A0
```

One-triple explicit-kernel smoke test:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/export_semialgebraic_system.py \
  --mode kernel --active-indices 0 \
  --out-prefix code/sage/out/smoke_kernel_A0
```

Export the known 10-active base extremal determinant system:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/export_semialgebraic_system.py \
  --mode det --known-label 'P(S(e,e,e),e,e,e)' --invert-d \
  --out-prefix code/sage/out/known_base_det_invd
```

Probe a determinant ideal with Singular, optionally over a finite field and with
generic affine linear sections:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/probe_determinant_ideals.py \
  --known-label 'P(S(e,e,e),e,e,e)' --invert-d \
  --characteristic 32003 --linear-sections 3 \
  --out-prefix code/sage/out/probe_known_base_invd_p32003_sec3
```

Compute a lex basis for a zero-dimensional exact slice and classify its real
roots:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/lexify_determinant_slice.py \
  --known-label 'P(S(e,e,e),e,e,e)' --invert-d --method modslimgb \
  --characteristic 0 --linear-sections 3 --linear-section-mode zero-sum-z \
  --seed 401 --out-prefix code/sage/out/lex_known_base_invd_QQ_zerosum401_sec3_z0

~/miniforge3/bin/mamba run -n sage python code/sage/classify_lex_slice_roots.py \
  --lex-json code/sage/out/lex_known_base_invd_QQ_zerosum401_sec3_z0.json
```

For a candidate exact slice certificate, check the sign of a principal minor at
algebraic real roots of the non-known factor:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/certify_lex_slice_minor.py \
  --lex-json code/sage/out/lex_known_base_invd_QQ_zerosum401_sec3_z0.json \
  --active-index 1 --minor-rows 0,2
```

To screen all active `2x2` principal minors on a lex slice:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/screen_lex_slice_minors.py \
  --lex-json code/sage/out/lex_known_base_invd_QQ_zerosum401_sec3_z0.json
```

To test whether the candidate relation
`minor((0,1,3); rows 0,2) + 46656*z0^2 = 0` selects the residual component in a
modular section probe:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/probe_minor_relation.py \
  --known-label 'P(S(e,e,e),e,e,e)' --invert-d \
  --characteristic 32003 --linear-sections 3 --linear-section-mode zero-sum-z \
  --seed 401 --out-prefix code/sage/out/probe_minor_relation_known_base_invd_p32003_zerosum401_sec3
```

To run the component-level separator/open-set tests:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/certify_component_relation.py \
  --known-label 'P(S(e,e,e),e,e,e)' --invert-d \
  --characteristic 32003 --out-prefix code/sage/out/cert_component_relation_known_base_invd_p32003_sep_z1mz0
```

Screen many subsets from a known active set:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/screen_determinant_subsets.py \
  --known-label 'P(S(e,e,e),e,e,e)' --subset-size 4 --max-subsets 0 \
  --out-prefix code/sage/out/screen_known_base_size4_all_p32003
```

Numerically sample real determinant sections through a known extremal and
classify the sampled roots by the active PSD and inactive inequalities:

```bash
.venv/bin/python code/sage/sample_real_determinant_sections.py \
  --known-label 'P(S(e,e,e),e,e,e)' --sections 3 --seed 101 --starts 800 \
  --out code/sage/out/real_section_known_base_seed101.json
```

The same sampler can center sections at a saved low-active margin point from
`verify/v30_margin_infimum.py`:

```bash
.venv/bin/python code/sage/sample_real_determinant_sections.py \
  --v30-input verify/out/v30_margin_6_78593.json \
  --sections 3 --section-mode random --seed 801 --starts 350 \
  --out code/sage/out/real_section_low_s6_78593_v30_seed801.json
```

This last script is a real/numerical probe, not a proof.  It is meant to guide
which semialgebraic constraints are worth certifying next.

The structured size-six low-active ansatz can be probed directly:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/probe_s6_ansatz.py \
  --characteristic 32003 --branch cplusd_sq5 --order lex \
  --out code/sage/out/s6_ansatz_cplusd_lex_p32003.json
```

After collecting several prime runs, reconstruct the residual q-factor by first
dividing the visible nuisance factor:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/reconstruct_univariate_crt.py \
  --variable q \
  --divide-factor '(q-3)*(q-5)^2*(q^2-3*q+6)^2*(q^2-9*q+24)^2' \
  --out code/sage/out/s6_ansatz_cplusd_q_residual_crt_6primes.json \
  code/sage/out/s6_ansatz_cplusd_lex_p32003.json \
  code/sage/out/s6_ansatz_cplusd_lex_p32009.json \
  code/sage/out/s6_ansatz_cplusd_lex_p32027.json \
  code/sage/out/s6_ansatz_cplusd_lex_p32029.json \
  code/sage/out/s6_ansatz_cplusd_lex_p32051.json \
  code/sage/out/s6_ansatz_cplusd_lex_p32057.json
```

Classify the resulting real q-values numerically:

```bash
.venv/bin/python code/sage/classify_s6_ansatz_numeric.py \
  --starts-per-q 800 \
  --out code/sage/out/classify_s6_ansatz_cplusd_numeric_s800.json
```

Exact Groebner/dimension probes should be run under an explicit timeout:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/export_semialgebraic_system.py \
  --mode det --known-label 'P(S(e,e,e),e,e,e)' --singular-compute \
  --out-prefix code/sage/out/known_base_det_compute

timeout 60s ~/miniforge3/bin/mamba run -n sage Singular -q \
  code/sage/out/known_base_det_compute.sing
```

The cofactor-patch model is mainly a direct transcription of the proof
condition.  For actual CAS attacks, the explicit-kernel model is currently the
better first target: it has more variables but much lower-degree equations.

See `RESULTS.md` for the local probe log, including the observed dimension
pattern for one-, two-, and three-active saturated determinant systems.
