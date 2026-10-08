# Tipping-point threshold for differential misclassification of the dichotomised functional outcome in neurocritical care trials

Analysis protocol, simulation code and complete output for a methodological study submitted to
*BMC Medical Research Methodology*.

**Registration.** OSF Registries, *Simulation Study Template* (ADEMP-PreReg):
[10.17605/OSF.IO/87ZJN](https://doi.org/10.17605/OSF.IO/87ZJN) — project at <https://osf.io/h8u35/>.

The registration is a **post-hoc deposit**: the analyses it describes were completed before it was
deposited, and it discloses that in its own text rather than leaving the date to be discovered.
`reports/estimator_eval_design.json` carries a timestamp that precedes the first output file of the
corresponding runs.

## What is here

| Path | Content |
|---|---|
| `scripts/simulate_tipping_point.py` | The closed form `Δ_tip = (se − f)(p₀ − p₁)/(1 − p₁)` and the layer-B boundary-error model |
| `scripts/calibrate_estimator_eval.py` | Type I error and detection rate over the 800-cell grid |
| `scripts/convergence_evidence.py` | Replicate doubling, grid refinement and `Δ_sig` resolution |
| `scripts/verify_bmcmrm_numbers.py` | Independent verifier: recomputes every reported number from the frozen JSON |
| `reports/estimator_eval_design.json` | The frozen analysis protocol — factors, values, legality assertions, assigned meaning of each output |
| `reports/estimator_eval_calibration.json` | Per-cell type I error and detection records |
| `reports/estimator_eval_delta_sig_coverage.json` | Per-cell `Δ_sig` and coverage records |
| `reports/estimator_eval_null_contour.json` | Layer-B contour records |
| `reports/estimator_eval_dgp_checks.json` | Data-generating-mechanism certification |
| `reports/estimator_eval_convergence.json` | The convergence block |
| `reports/estimator_eval_baseline.json` | Frozen baseline |
| `reports/simulation_tipping_point.json` | The earlier closed-form run |
| `reports/verify_bmcmrm_numbers_report.txt` | The verifier's clean-run report |

`MANIFEST.md` lists the byte length and MD5 of every file.

## Requirements

CPython 3.x, **standard library only** — `json`, `math`, `statistics`, `hashlib`, `random`, `re`,
`pathlib`, `itertools`, `argparse`. No third-party packages, no R, no compiled extensions, no network
access. Total runtime is minutes, not hours.

## Running the pipeline

Paths are relative to this directory, so run everything from the repository root:

```
python scripts/simulate_tipping_point.py     # -> reports/simulation_tipping_point.json
python scripts/calibrate_estimator_eval.py   # -> reports/estimator_eval_calibration.json
python scripts/convergence_evidence.py       # -> reports/estimator_eval_convergence.json
```

Reproducibility does not rest on the seeds alone. Cell streams are derived as a SHA-256 digest over a
namespaced cell key (for example `A|p0=0.45|or_true=1.0|se=0.7|f=0.0`) rather than consumed
sequentially from one generator, so inserting or reordering a cell cannot move any other cell's draws.
The master seed is `20261006`.

## Reproducing every reported number

```
python scripts/verify_bmcmrm_numbers.py            # 43 assertions over 31 recomputed quantities
python scripts/verify_bmcmrm_numbers.py --list     # the guard list and the mutation battery
python scripts/verify_bmcmrm_numbers.py --selftest # run all 28 deliberately corrupted inputs
```

This verifier recomputes each quantity directly from the frozen JSON rather than reading it out of the
manuscript, and separately checks that the manuscript *states* it. It therefore reads the manuscript
sources, which sit at `submission/bmcmrm_v1/` in the parent project and are not redistributed here
while the manuscript is under review. Place them at that path to reproduce the full run. The
`--selftest` battery is routed away from the delivery path so it cannot overwrite the shipped report,
and `VDR_OUTDIR` redirects the report elsewhere if you want to keep the frozen one intact:

```
VDR_OUTDIR=/tmp/vdr python scripts/verify_bmcmrm_numbers.py
```

## Scope and limits

- The study evaluates a **threshold**, not a new estimator. The comparison a reader would form from a
  finite trial is the object of study; the threshold itself is a closed-form quantity that requires no
  estimation from the data it is applied to.
- The threshold is **not recoverable from published arm risks**: those are *observed* risks, which
  already carry the misclassification in question. This is written into the protocol as a hard rule.
- **Convergence is not validity.** Stability under more repetitions, a finer grid and a resolved
  bisection bracket is reported, and a numerically stable estimator can still be the wrong estimator.
- No patient-level data are used. The inputs are the simulation parameters and the published search
  records reported in the manuscript's supplement.

## Licence

Code: MIT (see `LICENSE`). The registration text and the reported records are released under
CC-BY 4.0 through the OSF deposit.
