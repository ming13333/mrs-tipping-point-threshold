# MANIFEST

Every file below was deposited to OSF Registries on **2026-10-08** as registration
**10.17605/OSF.IO/87ZJN** and is reproduced here byte-for-byte. Total **8,344,154 bytes**.

| File | Bytes | MD5 |
|---|---:|---|
| `scripts/simulate_tipping_point.py` | 17,258 | `24ff3cbd407e2650a2b62498af36506e` |
| `scripts/calibrate_estimator_eval.py` | 34,430 | `de373eaf5590249c23cd81fd86d625fe` |
| `scripts/convergence_evidence.py` | 51,929 | `97b7bb127910f09ec8fb053be0057b15` |
| `scripts/verify_bmcmrm_numbers.py` | 73,289 | `35a1771146e182af273c7900594943e2` |
| `reports/estimator_eval_design.json` | 8,811 | `e4c37a9af26b262d6c6aef772f496ebc` |
| `reports/estimator_eval_baseline.json` | 15,069 | `729560bd5b8f77c31d006b95d9f5e2f4` |
| `reports/estimator_eval_calibration.json` | 3,648,054 | `9b7fac45536b0c724398ed4563a42f10` |
| `reports/estimator_eval_delta_sig_coverage.json` | 3,524,873 | `5d5542c940fb6418635e0b5bfaa2cb18` |
| `reports/estimator_eval_null_contour.json` | 127,933 | `ad2f26517db4945bd982e635ac950a70` |
| `reports/estimator_eval_dgp_checks.json` | 576,353 | `faccf6ba2f2793b2037f8abc8e2e6e5a` |
| `reports/estimator_eval_convergence.json` | 223,149 | `66d43a8cb8a58c2f2bb9c5945b9d4547` |
| `reports/simulation_tipping_point.json` | 36,755 | `f56e9c55f82388df3a3f2c0281b2d3b7` |
| `reports/verify_bmcmrm_numbers_report.txt` | 6,251 | `51596a0acb0033bf21b62c9b7fcdaef6` |

```bash
# verify the deposit is intact
md5sum -c <(awk -F'|' 'NR>5 && NF>3 {gsub(/[ `]/,"",$2); gsub(/[ `]/,"",$3); print $3"  "$2}' MANIFEST.md)
```

The registration's own copy of these files is at <https://osf.io/87zjn/files/>.
