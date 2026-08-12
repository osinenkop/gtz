# GTZ MLCore Watch Summary

- Updated: `2026-08-10T12:58:16+00:00`
- Started: `2026-08-10T08:49:46+00:00`
- Poll count: `26`
- Finished all jobs: `True`
- Deadline reached: `False`

## Job Status

| Job | State | Phase | Termination | Log signal |
| --- | --- | --- | --- | --- |
| `gtz-d6sup2-s7-p32029-c19-p3ih3w` | `SUCCEEDED` | `ENDED` | `JobCompleted` | defect 42; passing 153; `out/local_separator_s7_78612_omit_active1_D6_S6_32029_support19.json` |
| `gtz-d6sup2-s7-p32051-c19-9taq9k` | `SUCCEEDED` | `ENDED` | `JobCompleted` | defect 42; passing 153; `out/local_separator_s7_78612_omit_active1_D6_S6_32051_support19.json` |
| `gtz-d6sup2-s7-p32057-c19-6njtnj` | `SUCCEEDED` | `ENDED` | `JobCompleted` | defect 42; passing 153; `out/local_separator_s7_78612_omit_active1_D6_S6_32057_support19.json` |
| `gtz-d6sup2-s8-p32029-c6-oxrgpb` | `SUCCEEDED` | `ENDED` | `JobCompleted` | defect 63; passing 262; `out/local_separator_s8_79656_omit_active0_D6_S6_32029_support6.json` |
| `gtz-d6sup2-s8-p32051-c6-eb1njn` | `SUCCEEDED` | `ENDED` | `JobCompleted` | defect 63; passing 262; `out/local_separator_s8_79656_omit_active0_D6_S6_32051_support6.json` |
| `gtz-d6sup2-s8-p32057-c6-h11gtp` | `SUCCEEDED` | `ENDED` | `JobCompleted` | defect 63; passing 262; `out/local_separator_s8_79656_omit_active0_D6_S6_32057_support6.json` |

## Separator JSONs

| File | p | D | S | rank | gen rank | defect | kernel | rank-only | candidates | passing | max |s(root)| |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | --- |
| `local_separator_s7_78612_omit_active1_D5_S4_p32003.json` | 32003 | 5 | 4 | 30916 | 29920 | 5 | 115 | False | 12 | 12 |  |
| `local_separator_s7_78612_omit_active1_D6_S6_32003.json` | 32003 | 6 | 6 | 87073 | 79399 | 334 | 1015 | False | 834 | 832 |  |
| `local_separator_s7_78612_omit_active1_D6_S6_32029_support19.json` | 32029 | 6 | 6 | 79629 | 79399 | 42 | 723 | False | 153 | 153 | 243681.74059427327154054882885177276660175945903724472285754 |
| `local_separator_s8_79656_omit_active0_D5_S4_p32003.json` | 32003 | 5 | 4 | 30860 | 29865 | 6 | 171 | False | 21 | 21 |  |
| `local_separator_s8_79656_omit_active0_D6_S6_32003.json` | 32003 | 6 | 6 | 86680 | 79071 | 399 | 1408 | False | 1408 | 1408 |  |
| `local_separator_s7_78612_omit_active1_D6_S6_32051_support19.json` | 32051 | 6 | 6 | 79629 | 79399 | 42 | 723 | False | 153 | 153 | 287169.11908887968016568133169204723431392072423823416247836 |
| `local_separator_s7_78612_omit_active1_D6_S6_32057_support19.json` | 32057 | 6 | 6 | 79629 | 79399 | 42 | 723 | False | 153 | 153 | 204726.64892714707308149940009404175371448776025083175089560 |
| `local_separator_s8_79656_omit_active0_D6_S6_32029_support6.json` | 32029 | 6 | 6 | 79327 | 79071 | 63 | 1072 | False | 262 | 262 | 327895.13343823753544397828221745281897762478536019086817632 |
| `local_separator_s8_79656_omit_active0_D6_S6_32051_support6.json` | 32051 | 6 | 6 | 79327 | 79071 | 63 | 1072 | False | 262 | 262 | 420251.52224231559135110700452412759039667548732066145357773 |
| `local_separator_s8_79656_omit_active0_D6_S6_32057_support6.json` | 32057 | 6 | 6 | 79327 | 79071 | 63 | 1072 | False | 262 | 262 | 344355.64788303101827919285902621653761932073275332404968368 |

## Conclusion

Potential breakthrough: at least one full extraction has a nonzero-at-root separator candidate.
