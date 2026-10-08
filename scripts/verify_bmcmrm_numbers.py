# -*- coding: utf-8 -*-
"""verify_bmcmrm_numbers.py -- bmcmrm_v1 reconciliation.

Reconciles every number that appears in the manuscript, the eMethods supplement
and the eTable 1 supplement against the artefacts that produced it.

Sources of truth, in order of precedence:

  1. `scripts/simulate_tipping_point.py` -- the legacy closed form and the
     boundary-error model. Imported, never retyped, so the manuscript and the
     code cannot drift apart silently.
  2. `reports/simulation_tipping_point.json` -- the legacy run.
  3. The six W1 artefacts in `reports/estimator_eval_*.json`.

Three families of check:

  A. every numeric token in prose must be *accounted for* -- either structural
     (a citation number, a section number, a year), an external figure carried
     over from a cited publication, or a value this repository derives. A token
     that belongs to none of the three fails the run. This is a whitelist, so a
     number added to the manuscript without a decision is caught.
  B. every derived value must be *recomputed* and match the text to the decimal
     place the manuscript prints.
  C. the W1 findings that v12 does not yet contain must be registered as
     pending, and v12 must not assert anything stronger than v12 can support.

Discipline notes carried over from steps 1-6:

  * forbidden-pattern guards are built from fragments so they cannot match
    themselves;
  * each guard must be *reachable*: `--selftest` injects one mutation per guard
    and requires the run to fail, so a guard that cannot fail is reported as
    inert rather than as a pass;
  * the tolerance used for width and residual comparisons is the double
    precision floor, 2**-52, never an ideal value derived from an iteration
    count.

Usage:
    python scripts/verify_manuscript_numbers.py            # the real run
    python scripts/verify_manuscript_numbers.py --list     # guards + injections
    python scripts/verify_manuscript_numbers.py --selftest # mutation battery

Exit code 0 only if every guard passes AND, under --selftest, every injection
is caught.
"""
from __future__ import annotations

import io
import json
import os
import re
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import simulate_tipping_point as legacy  # noqa: E402  the closed form itself

SUBMISSION = ROOT / "submission" / "bmcmrm_v1"
MS = SUBMISSION / "manuscript_submission.md"
EM = SUBMISSION / "esupplement_eMethods.md"
ET = SUBMISSION / "esupplement_eTable1_search_strategy.md"
AB = SUBMISSION / "abstract.md"
CP = SUBMISSION / "reporting_compliance.md"
LEGACY_JSON = ROOT / "reports" / "simulation_tipping_point.json"

# The journal's own limit for a structured abstract at BMC Medical Research
# Methodology.  It is a hard submission gate, so it is asserted, not eyeballed.
ABSTRACT_LIMIT = 350

SPEC = ROOT / "reports" / "estimator_eval_design.json"
BASE = ROOT / "reports" / "estimator_eval_baseline.json"
DGP = ROOT / "reports" / "estimator_eval_dgp_checks.json"
CAL = ROOT / "reports" / "estimator_eval_calibration.json"
DSIG = ROOT / "reports" / "estimator_eval_delta_sig_coverage.json"
CONT = ROOT / "reports" / "estimator_eval_null_contour.json"
CONV = ROOT / "reports" / "estimator_eval_convergence.json"

# Injection runs must never write into the delivery tree: the self-test spawns
# one subprocess per mutation and the last one to finish wins, so a hardcoded
# path leaves the shipped report recording a *failing* run.  The output
# directory is therefore routable, and the self-test routes it to scratch.
OUTDIR = Path(os.environ.get("VDR_OUTDIR") or (ROOT / "reports"))
OUTDIR.mkdir(parents=True, exist_ok=True)
REPORT = OUTDIR / "verify_bmcmrm_numbers_report.txt"

MACHINE_FLOOR = 2.0 ** -52
TOL = 1e-9

INJECT = os.environ.get("VMN_INJECT", "")

# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------
_lines: list[str] = []
_guards: list[tuple[str, bool, str]] = []


def guard(name: str, ok: bool, detail: str = "") -> bool:
    _guards.append((name, bool(ok), detail))
    line = f"[{'PASS' if ok else 'FAIL'}] {name}"
    if detail:
        line += f"  ->  {detail}"
    _lines.append(line)
    print(line)
    return bool(ok)


def note(text: str) -> None:
    _lines.append(f"       {text}")
    print(f"       {text}")


# ---------------------------------------------------------------------------
# injection table: one mutation per guard family
# ---------------------------------------------------------------------------
INJECTIONS: dict[str, tuple[str, str, str, str]] = {
    # (file key, old, new, human description)
    "ms_table1_cell": ("ms", "| 15.4      | 17.6      | 19.8      | 21.6      |",
                       "| 15.4      | 17.6      | 19.9      | 21.6      |",
                       "corrupt one Table 1 cell"),
    "ms_table1_p1_label": ("ms", "0.60 (p₁ = 0.423)", "0.60 (p₁ = 0.424)",
                           "corrupt one p1 annotation in Table 1"),
    "ms_table1_drop_cell": ("ms", "| 0.60              | 0.80 (p₁ = 0.545)    | 8.4       | 9.6       | 10.8      | 11.8      |",
                            "",
                            "delete one Table 1 row"),
    "ms_headline": ("ms", "Δ_tip = 19.8 percentage points", "Δ_tip = 19.9 percentage points",
                    "corrupt the headline threshold in the abstract"),
    "ms_se_scan": ("ms", "15.4–21.6", "15.4–21.7", "corrupt the se scan range"),
    "ms_two_way": ("ms", "spanning 17.7 to 26.0 percentage points",
                   "spanning 17.7 to 27.0 percentage points", "corrupt the two-way family range"),
    "ms_upgrade_rates": ("ms", "upgrade rate of 71% in the index arm against 15% in the reference arm",
                         "upgrade rate of 61% in the index arm against 15% in the reference arm",
                         "corrupt the index-arm upgrade rate"),
    "ms_artefact": ("ms", "the apparent asymmetry runs from nothing to about 4.5 percentage points",
                    "the apparent asymmetry runs from nothing to about 4.6 percentage points",
                    "corrupt the common-process artefact in the abstract"),
    "ms_substitution": ("ms", "obtains 16.6 percentage points where the correct figure is 19.8, an underestimate of 16%",
                        "obtains 16.6 percentage points where the correct figure is 19.8, an underestimate of 6%",
                        "corrupt the observed-substitution underestimate"),
    "ms_observed_risks": ("ms", "the observed risks are 0.495 and 0.381",
                          "the observed risks are 0.495 and 0.381", "placeholder, replaced below"),
    "em_grid_size": ("em", "(16 levels) and b₁ from 0 to 0.40 in steps of 0.05 (9 levels), giving 144 combinations",
                   "(16 levels) and b₁ from 0 to 0.40 in steps of 0.05 (9 levels), giving 72 combinations",
                   "reintroduce the grid-size error this step was written to catch"),
    "em_mc_cell": ("em", "| **1.0** | **19.80** | **1.0000** | **1.0000** | 0.841 – 1.197 |",
                   "| **1.0** | **19.80** | **1.0001** | **1.0000** | 0.841 – 1.197 |",
                   "corrupt the Monte-Carlo median at Delta_tip"),
    "em_contour_cell": ("em", "| 0.25 | 0.7949 | 22.90 | 0.9291 |",
                        "| 0.25 | 0.7949 | 22.91 | 0.9291 |",
                        "corrupt one OR=1 contour row"),
    "em_part2_cell": ("em", "| Reference | 0.15 | 0.10 | 0.9709 | 0.0700 | 0.5655 |",
                      "| Reference | 0.15 | 0.10 | 0.9709 | 0.0701 | 0.5655 |",
                      "corrupt the reference-arm effective false-positive rate"),
    "em_artefact_max": ("em", "| \\|apparent excess\\|, maximum | **4.50 pp** |",
                        "| \\|apparent excess\\|, maximum | **4.60 pp** |",
                        "corrupt the artefact maximum"),
    "em_ratio": ("em", "| Ratio of the two | **0.997** |",
                 "| Ratio of the two | **0.999** |", "corrupt the closed-form ratio"),
    "em_slice_endpoint": ("em", "from 0.5606 at a₁ = 0.05 to 1.0813 at a₁ = 0.80",
                          "from 0.5606 at a₁ = 0.05 to 1.0814 at a₁ = 0.80",
                          "corrupt the one-way slice upper endpoint"),
    "etable1_yield": ("et", "| `(\"modified Rankin Scale\"[tiab] OR mRS[tiab]) AND (misclassification[tiab] OR misclassified[tiab] OR \"measurement error\"[tiab] OR interobserver[tiab] OR \"inter-rater\"[tiab] OR reliability[tiab] OR agreement[tiab])` | 1,188 |",
                      "| `(\"modified Rankin Scale\"[tiab] OR mRS[tiab]) AND (misclassification[tiab] OR misclassified[tiab] OR \"measurement error\"[tiab] OR interobserver[tiab] OR \"inter-rater\"[tiab] OR reliability[tiab] OR agreement[tiab])` | 1,290 |",
                      "make the eTable 1 S1 yield disagree with the manuscript"),
    "unregistered_number": ("ms", "It fixes the quantity to be reported",
                            "It fixes the quantity to be reported, which came to 17.3 in our data",
                            "insert a number nobody accounted for"),
    "overclaim_ci": ("ms", "The 19.8 figure is the point at which the point estimate is null",
                     "The 19.8 figure is the point at which the confidence interval becomes null",
                     "assert that Delta_tip is a significance threshold"),
    "conflate_failure": ("ms", "the asymmetry is not free to control",
                         "the closed form is wrong by up to 4.5 percentage points, and the asymmetry is not free to control",
                         "conflate the apparent asymmetry with the closed-form failure"),
    "smallest_claim": ("ms", "it is the slice obtained when the downgrade rate is held equal",
                       "it is the smallest value in the family, obtained when the downgrade rate is held equal",
                       "claim the one-way value is the family minimum"),
    "significance_claim": ("ms", "the interval around it is a separate quantity",
                           "the interval around it is a separate quantity at the 5 per cent level, where significance is reached",
                           "assert a significance threshold"),
    "ms_type1_swapped": ("ms", "0.0525", "0.0492",
                         "swap the pooled type I error for another whitelisted value"),
    "one_in_five": ("ms", "Two features of Table 1 matter.",
                    "Roughly one in five outcome assignments would have to be wrong. Two features of Table 1 matter.",
                    "reintroduce the share-of-assignments misstatement"),
}
INJECTIONS["ms_observed_risks"] = (
    "ms", "the observed risks are 0.495 and 0.381", "the observed risks are 0.496 and 0.381",
    "corrupt the observed reference risk")
INJECTIONS["cp_wordcount_stale"] = (
    "cp", "**350 words**", "**351 words**",
    "leave one of the two abstract word-count declarations behind")
INJECTIONS["ab_abstract_diverges"] = (
    "ab", "Two threats to the dichotomised functional outcome in neurocritical care trials are conflated",
    "Two threats to the functional outcome are conflated",
    "let the standalone abstract drift away from the manuscript copy")
INJECTIONS["ms_abstract_over_limit"] = (
    "ms", "does not replace blinded ascertainment; what the audit establishes is that the comparison it would be read against is essentially never reported.",
    "does not replace blinded ascertainment and is not a substitute for the prospective measurement of agreement between the two arms; what the audit establishes is that the comparison it would be read against is essentially never reported in any of the records we screened.",
    "push the abstract past the journal's 350-word gate")


def load_texts() -> dict[str, str]:
    t = {"ms": MS.read_text(encoding="utf-8"),
         "em": EM.read_text(encoding="utf-8"),
         "et": ET.read_text(encoding="utf-8"),
         "ab": AB.read_text(encoding="utf-8"),
         "cp": CP.read_text(encoding="utf-8")}
    if INJECT:
        key, old, new, _desc = INJECTIONS[INJECT]
        if old not in t[key]:
            raise SystemExit(f"INJECTION ERROR: anchor for {INJECT!r} absent from {key}")
        t[key] = t[key].replace(old, new, 1)
    return t


# ---------------------------------------------------------------------------
# text plumbing
# ---------------------------------------------------------------------------
CODE_FENCE = re.compile(r"```.*?```", re.S)
CITE = re.compile(r"\[[0-9,\s–\-]+\]")
LABEL = re.compile(r"\*\*(?:Table \d|Fig\. \d|Box \d)[^*]*\*\*")
SEP_ROW = re.compile(r"^\|[\s|\-:]*\|$", re.M)
YEAR = re.compile(r"\b(?:19|20)\d\d\b")
ORCID = re.compile(r"ORCID:? \S+")
POSTCODE = re.compile(r"Lishui \d+")
TOOLVER = re.compile(r"WorkBuddy version [\d.]+, powered by the [\w.\-]+ model")
SEED = re.compile(r"\b20\d{6}\b")
PYVER = re.compile(r"Python 3\.\d+")
MONEY = re.compile(r"£[\d,]+")
COMMA_NUM = re.compile(r"\b\d{1,3},\d{3}\b")
TOKEN = re.compile(r"(?<![\w.])\d+(?:\.\d+)?(?![\w.])")


def strip_noise(text: str, drop_after: str | None = None) -> str:
    if drop_after:
        text = text[: text.index(drop_after)]
    text = CODE_FENCE.sub("", text)
    text = CITE.sub("", text)
    text = LABEL.sub("", text)
    text = SEP_ROW.sub("", text)
    text = ORCID.sub("", text)
    text = POSTCODE.sub("", text)
    text = TOOLVER.sub("", text)
    text = SEED.sub("", text)
    text = PYVER.sub("", text)
    text = MONEY.sub("", text)
    text = COMMA_NUM.sub("", text)
    text = YEAR.sub("", text)
    return text


def d(value: float, places: int) -> str:
    """Round half up at `places`, the way a typesetter would print it."""
    quant = Decimal(1).scaleb(-places)
    return str(Decimal(repr(float(value))).quantize(quant, rounding=ROUND_HALF_UP))


def tokens(text: str) -> dict[str, set[int]]:
    out: dict[str, set[int]] = {}
    for i, line in enumerate(strip_noise(text).splitlines(), 1):
        for m in TOKEN.finditer(line):
            out.setdefault(m.group(0), set()).add(i)
    return out


# ---------------------------------------------------------------------------
# A. accounting tables
# ---------------------------------------------------------------------------
# Structural: not scientific quantities. Ordinals, citation counts, section and
# table numbers, dates, ages, and the small integers that carry grammar.
STRUCTURAL_MS = {
    "0": "the f = 0 anchor; the constraint 0 <= f + Delta <= 1; zero events",
    "1": "one-way; 'one line'; OR = 1; Table 1; Box 1; single direction",
    "2": "the mRS 2/3 boundary; two arms; two directions of error",
    "3": "the mRS 3-6 dichotomy; two directions",
    "4": "four decades of reliability work; four sensitivities in Table 1",
    "5": "five per cent alpha; five decades",
    "6": "the mRS 3-6 dichotomy; mRS 6 in the abbreviations",
    "7": "abbreviation list ordinals",
    "8": "eight grade levels in the eMethods distributions",
    "9": "no bare 9 survives comma stripping",
    "15": "the reference-arm upgrade rate is derived; see DERIVED_MS",
    "16": "the 16 per cent underestimate is derived; see DERIVED_MS",
    "17": "17 shape combinations (external, ref 17 family) and the two-way lower end",
    "20": "twenty-two studies split; see 22",
    "22": "22 studies of the stroke severity scale (external, ref 27)",
    "24": "the search date, 24 September 2026",
    "26": "26 intensive care units (external, ref 1)",
    "46": "a meta-analysis of 46 studies (external, ref 19)",
    "60": "agreement in about 60 per cent of patients (external, ref 18)",
    "71": "the index-arm upgrade rate is derived; see DERIVED_MS",
    "90": "the 90-day endpoint; se = 0.90 is a parameter",
    "95": "95 per cent confidence interval; 95 per cent MC interval",
    "100": "the bound f + Delta <= 1 written as 100 x in Table 1's note; GBP 100,000",
    "189": "1,189 records; asserted equal to eTable 1 S1",
    "205": "9,205 patients (external, ref 29)",
    "295": "GBP 2,295 per event (external, ref 36)",
    "000": "GBP 100,000 (external, ref 36)",
    "704": "704 ventilated adults (external, ref 1)",
    "860": "an 860-patient TBI trial dataset (external, ref 28)",
    "179": "ORCID", "0005": "ORCID", "0009": "ORCID", "4834": "ORCID",
    "323000": "postcode",
    "1980": "since the 1980s",
    "2026": "the search year",
}

# External figures carried over from cited publications. Each is an input to
# the manuscript's argument, not an output of this repository.
EXTERNAL_MS = {
    "10.6": "share of scores classified incorrectly, van der Ende 2022 (ref 16)",
    "25.8": "local-central disagreement at 90 days, Sluis 2025 (ref 14)",
    "50.3": "poor outcome at 90 days, MAC-ICU (ref 1)",
    "0.97": "adjusted OR, MAC-ICU (ref 1)",
    "0.48": "CI lower bound, MAC-ICU (ref 1)",
    "1.96": "CI upper bound, MAC-ICU (ref 1)",
    "0.68": "kappa, Sluis 2025 (ref 14); also the upper end of the unweighted range in ref 17",
    "0.87": "weighted kappa, Shawkat 2025 (ref 15)",
    "0.56": "weighted kappa lower end, Nair 2025 SR (ref 17)",
    "0.27": "unweighted kappa lower end, Nair 2025 SR (ref 17)",
    "0.65": "pooled inter-rater reliability, Tvrda 2025 meta-analysis (ref 19)",
    "2.1": "per cent of participants misclassified differentially in randomised stroke trials (ref 30)",
    "6": "the upper end of that same 2.1 to 6 per cent range (ref 30)",
    "9": "9,205-patient dataset; the 9 is absorbed by comma stripping",
}

# Scenario parameters fixed by the manuscript.
PARAM_MS = {
    "0.50": "reference-risk scenario",
    "0.55": "reference-risk scenario (headline)",
    "0.60": "hypothesised true OR scenario",
    "0.70": "sensitivity scenario / hypothesised OR scenario",
    "0.80": "sensitivity scenario / hypothesised OR scenario",
    "0.90": "headline sensitivity",
    "0.98": "upper sensitivity scenario",
    "0.05": "0.05 step in the reference-risk sensitivity statement",
}

# Table 1 cells and every other derived value printed in the manuscript.
DERIVED_MS = {
    "0.375": "p1 implied by OR 0.60 at p0 0.50",
    "0.412": "p1 implied by OR 0.70 at p0 0.50",
    "0.423": "p1 implied by OR 0.60 at p0 0.55",
    "0.444": "p1 implied by OR 0.80 at p0 0.50",
    "0.461": "p1 implied by OR 0.70 at p0 0.55",
    "0.474": "p1 implied by OR 0.60 at p0 0.60",
    "0.494": "p1 implied by OR 0.80 at p0 0.55",
    "0.512": "p1 implied by OR 0.70 at p0 0.60",
    "0.545": "p1 implied by OR 0.80 at p0 0.60",
    "0.381": "observed index risk under the headline scenario",
    "0.495": "observed reference risk under the headline scenario",
    "14.0": "Table 1", "16.0": "Table 1", "18.0": "Table 1", "19.6": "Table 1",
    "10.5": "Table 1", "12.0": "Table 1", "13.5": "Table 1", "14.7": "Table 1",
    "7.0": "Table 1", "8.0": "Table 1", "9.0": "Table 1", "9.8": "Table 1",
    "15.4": "Table 1", "17.6": "Table 1", "19.8": "headline threshold and Table 1",
    "21.6": "Table 1 and the se scan",
    "11.6": "Table 1", "13.2": "Table 1", "14.9": "Table 1", "16.2": "Table 1",
    "7.7": "Table 1", "8.8": "Table 1", "9.9": "Table 1", "10.8": "Table 1",
    "16.8": "Table 1", "19.2": "Table 1", "23.5": "Table 1",
    "12.6": "Table 1", "14.4": "Table 1",
    "8.4": "Table 1", "9.6": "Table 1", "11.8": "Table 1",
    "4.5": "common-process apparent-asymmetry maximum",
    "16.6": "threshold computed on observed rather than true risks",
    "8.5": "same-denominator conversion of the 10.6% figure",
    "10.1": "same-denominator conversion of the 10.6% figure",
    "20.6": "same-denominator conversion of the 25.8% figure",
    "24.5": "same-denominator conversion of the 25.8% figure",
    "71": "index-arm upgrade rate for a null",
    "26.0": "two-way family upper end, printed to one decimal",
    "17.7": "two-way family lower end",
    # --- sampling-behaviour block (Part 5/6 of the eMethods) ---
    "160": "calibration cells and replication-doubling cells",
    "640": "Delta_sig cells",
    "800": "coverage cells",
    "68": "layer-B combinations with an OR = 1 contour",
    "73": "layer-B combinations in total",
    "200": "coarse-grid cell count",
    "630": "refined-grid cell count",
    "2.13": "lower end of the signed layer-B closed-form failure",
    "6.19": "upper end of the signed layer-B closed-form failure",
    "1.91": "mean Delta_sig minus Delta_tip, percentage points",
    "0.0525": "pooled type I error at the true null",
    "0.0008": "Monte Carlo standard error of that pooled rate",
    "0.0492": "lower end of the type I error over the p0 x se plane",
    "0.0544": "upper end of the type I error over the p0 x se plane",
    "0.722": "mean Delta_sig/Delta_tip at 100 per arm",
    "0.934": "mean Delta_sig/Delta_tip at 2,000 per arm",
    "87.8": "pooled uncorrected coverage",
    "95.2": "pooled corrected coverage",
    "94.2": "uncorrected coverage at 100 per arm",
    "79.6": "uncorrected coverage at 2,000 per arm",
    "0.023": "mean paired z on doubling replication",
    "1.005": "SD of that paired z",
    "0.124": "largest movement of Delta_sig on refined bisection",
    "0.00008": "movement of the pooled type I error on grid refinement",
    "2.8": "largest residual of the closed-form failure identity, x 1e-14",
    "10": "the exponent 10 in 2.8 x 10^-14",
    "430": "new cells added by midpoint refinement (630 - 200)",
}

# eMethods accounting.
STRUCTURAL_EM = {
    "0": "a and b anchor rows; zero events",
    "1": "one-way; P(observed >= 3); 2.5th and 97.5th percentiles",
    "2": "the mRS 2/3 boundary",
    "3": "Python 3.13; Part 3",
    "4": "Part 4",
    "5": "the 0.05 grid step",
    "6": "Part 6; mRS 6",
    "7": "reproduction section",
    "9": "no bare 9 survives comma stripping",
    "15": "the reference-arm upgrade rate is a parameter; see PARAM_EM",
    "16": "16 shape combinations",
    "32": "reference numbers in the closing paragraph",
    "34": "reference numbers in the closing paragraph",
    "35": "reference numbers in the closing paragraph (Pekkanen 2006)",
    "37": "reference numbers in the closing paragraph (Ma 2025)",
    "5.1": "section number", "5.2": "section number", "5.3": "section number",
    "5.4": "section number", "6.1": "section number", "6.2": "section number",
    "6.3": "section number", "6.4": "section number", "7.1": "section number",
    "7.2": "section number", "7.3": "section number", "7.4": "section number",
    "0.45": "lowest reference risk in the calibration grid",
    "0.65": "highest reference risk in the calibration grid",
    "0.70": "sensitivity level", "0.98": "sensitivity level",
    "0.75": "the 0.75 pp classification band",
    "1.03": "the 0.5 pp band expressed in paired standard errors",
    "0.0517": "pooled type I error at two decimals, Part 6",
    "0.0525": "pooled type I error at four decimals, Part 5",
    "100": "the factor 100 in the failure identity; 100 patients per arm",
    "256": "SHA-256",
    "95": "95 per cent MC interval, and the reference list of section 5",
    "5.3": "tool model name",
    "5.5": "tool version",
}
EXTERNAL_EM: dict[str, str] = {}
PARAM_EM = {
    "0.15": "reference-arm upgrade rate a0, and the common-process rate",
    "0.10": "reference-arm downgrade rate b0, and b1 held equal to it",
}
DERIVED_EM = {
    "0.90": "headline sensitivity",
    "0.55": "reference true risk",
    "0.60": "hypothesised true OR, and the upper concentration value",
    "0.423": "index true risk",
    "0.5998": "true OR implied by the two grade distributions",
    "0.090": "grade distribution", "0.150": "grade distribution",
    "0.210": "grade distribution", "0.160": "grade distribution",
    "0.110": "grade distribution", "0.130": "grade distribution",
    "0.200": "grade distribution", "0.217": "grade distribution",
    "0.120": "grade distribution", "0.093": "grade distribution",
    "0.05": "grid step", "0.80": "grid upper end and the one-way slice upper end",
    "0.40": "grid upper end", "0.35": "grid step",
    "0.30": "grid step and concentration value",
    "0.50": "concentration value",
    "0.00": "grid lower end",
    "0.20": "contour row label, b1 = 0.20",
    "0.25": "contour row label, b1 = 0.25",
    "19.80": "Delta_tip to two decimals",
    "9.90": "Monte-Carlo row", "15.84": "Monte-Carlo row", "23.76": "Monte-Carlo row",
    "29.70": "Monte-Carlo row",
    "0.7947": "analytic observed OR", "0.7977": "MC median",
    "0.668": "MC 2.5th percentile", "0.949": "MC 97.5th percentile",
    "0.9126": "analytic observed OR", "0.9120": "MC median",
    "0.767": "MC 2.5th percentile", "1.088": "MC 97.5th percentile",
    "1.0000": "analytic OR and MC median at Delta_tip; also index se in the contour table",
    "0.841": "MC 2.5th percentile", "1.197": "MC 97.5th percentile",
    "1.0957": "analytic observed OR", "1.0964": "MC median",
    "0.923": "MC 2.5th percentile", "1.292": "MC 97.5th percentile",
    "1.2576": "analytic observed OR", "1.2572": "MC median",
    "1.058": "MC 2.5th percentile", "1.491": "MC 97.5th percentile",
    "0.5": "Monte-Carlo row label, the multiple of Delta_tip",
    "0.8": "Monte-Carlo row label, the multiple of Delta_tip",
    "1.0": "Monte-Carlo row label, the multiple of Delta_tip",
    "1.2": "Monte-Carlo row label, the multiple of Delta_tip",
    "1.5": "Monte-Carlo row label, the multiple of Delta_tip",
    "0.9709": "reference effective se", "0.0700": "reference effective f",
    "0.5655": "reference observed risk",
    "0.9716": "index effective se", "0.0564": "index effective f",
    "0.4436": "index observed risk",
    "1.36": "difference between the two effective false-positive rates",
    "0.712": "index upgrade rate for a null",
    "71": "index upgrade rate as a percentage",
    "19.8": "headline threshold quoted in prose",
    "19.83": "closed-form Delta_tip at the reference effective (se, f)",
    "19.78": "effective excess actually delivered at the crossing",
    "0.997": "ratio of the two",
    "0.2678": "index effective false-positive rate at the crossing",
    "0.3": "the two agree to within 0.3 per cent",
    "0.5606": "observed OR at a1 = 0.05 on the one-way slice",
    "1.0813": "observed OR at a1 = 0.80 on the one-way slice",
    "0.6567": "contour a1 at b1 = 0.00", "17.70": "contour effective excess",
    "0.6843": "contour a1 at b1 = 0.05", "18.74": "contour effective excess",
    "0.9858": "contour index se at b1 = 0.05",
    "0.7120": "contour a1 at b1 = 0.10",
    "0.7396": "contour a1 at b1 = 0.15", "20.82": "contour effective excess",
    "0.9574": "contour index se at b1 = 0.15",
    "0.7673": "contour a1 at b1 = 0.20", "21.86": "contour effective excess",
    "0.9433": "contour index se at b1 = 0.20",
    "0.7949": "contour a1 at b1 = 0.25", "22.90": "contour effective excess",
    "0.9291": "contour index se at b1 = 0.25",
    "0.8226": "contour a1 at b1 = 0.30", "23.94": "contour effective excess",
    "0.9149": "contour index se at b1 = 0.30",
    "0.8502": "contour a1 at b1 = 0.35", "24.98": "contour effective excess",
    "0.9007": "contour index se at b1 = 0.35",
    "0.8779": "contour a1 at b1 = 0.40", "26.02": "contour effective excess",
    "0.8865": "contour index se at b1 = 0.40",
    "4.50": "apparent-excess maximum",
    "23.2": "maximum as a fraction of the threshold",
    "72": "retired: the supplement used to claim 72 combinations here",
    "16": "the a1 level count, and the 16 shape combinations",
    "9": "the b1 level count",
    "144": "the two-way grid size, 16 a1 levels by 9 b1 levels",
    # --- Part 5: sampling behaviour ---
    "160": "calibration cells; replication-doubling cells; the shared-cell count",
    "20": "summarising cells on the baseline-risk by sensitivity plane",
    "0.05252": "pooled type I error to five decimals",
    "0.00077": "Monte Carlo standard error of that pooled rate",
    "0.0370": "lowest single-cell type I error",
    "0.0665": "highest single-cell type I error",
    "0.0109": "largest single-cell Monte Carlo 95 per cent half-width",
    "0.0492": "plane lower end", "0.0544": "plane upper end",
    "0.0052": "spread of the type I error over the plane",
    "640": "Delta_sig cells",
    "8,000": "replicates per Delta_sig and coverage cell",
    "26": "bisection steps used to locate Delta_sig",
    "0.8583": "grid-wide mean Delta_sig/Delta_tip",
    "0.3535": "lowest Delta_sig/Delta_tip", "0.9727": "highest Delta_sig/Delta_tip",
    "1.91": "mean shortfall of Delta_sig below Delta_tip",
    "5.63": "largest shortfall of Delta_sig below Delta_tip",
    "0.7219": "mean ratio at 100 per arm", "0.8700": "mean ratio at 500 per arm",
    "0.9070": "mean ratio at 1,000 per arm", "0.9341": "mean ratio at 2,000 per arm",
    "0.7175": "lowest ratio at 500 per arm", "0.9465": "highest ratio at 500 per arm",
    "0.7978": "lowest ratio at 1,000 per arm", "0.9629": "highest ratio at 1,000 per arm",
    "0.8578": "lowest ratio at 2,000 per arm",
    "0.9063": "highest ratio at 100 per arm",
    "2.8936": "c in the 1 - c/sqrt(n) fit", "2.66": "coefficient of variation of that fit",
    "1.56": "largest residual of that fit, per cent",
    "800": "coverage cells", "0.8779": "pooled uncorrected coverage",
    "0.0011": "Monte Carlo half-width of that coverage",
    "0.9516": "pooled corrected coverage",
    "0.9420": "uncorrected coverage at 100 per arm", "0.9075": "at 500 per arm",
    "0.8661": "uncorrected coverage at 1,000 per arm", "0.7959": "at 2,000 per arm",
    "0.9513": "corrected coverage at 100 per arm", "0.9518": "corrected coverage at 500 and 1,000 per arm",
    "0.9515": "corrected coverage at 2,000 per arm",
    "73": "layer-B combinations in total", "68": "combinations with a contour",
    "1.02": "lower end of the inadmissible algebraic solution", "1.15": "upper end",
    "2.13": "lower end of the signed closed-form failure",
    "6.19": "upper end of the signed closed-form failure",
    "2.8": "largest residual of the failure identity, x 1e-14",
    "10": "the exponent 10 in 2.8 x 10^-14 and 1.6 x 10^-9",
    # --- Part 6: convergence ---
    "33": "cells whose value moved when an extra cell was inserted", "40": "cells in the demonstration block",
    "0.040": "value before the insertion", "0.048": "value after the insertion",
    "0.0230": "mean paired z on doubling", "1.0046": "SD of that paired z",
    "0.7979": "expected mean absolute z", "0.8230": "observed mean absolute z",
    "6": "cells beyond two paired standard errors", "7.3": "expected number beyond two",
    "107": "cells within one paired standard error", "154": "cells within two",
    "1.73": "sqrt(3), the unpaired-to-paired tolerance inflation",
    "39": "error of the unpaired prediction, per cent",
    "1.5": "error of the paired prediction, per cent",
    "0.003497": "empirical SD of 4,000 nested pairs", "0.003446": "paired prediction",
    "0.004873": "unpaired prediction",
    "1.0": "the classification band, in percentage points",
    "17": "cells whose classification changed on doubling",
    "48": "cells changing at a 0.5 pp band", "23": "cells changing at a 0.75 pp band",
    "2": "cells changing at a 1.5 pp band",
    "200": "coarse-grid cells", "630": "refined-grid cells",
    "5.85": "lower end of the Delta_tip range", "31.85": "upper end of the Delta_tip range",
    "0.051678": "pooled type I error, coarse grid",
    "0.051600": "pooled type I error, refined grid",
    "0.000078": "their difference", "0.001199": "combined Monte Carlo half-width",
    "0.011": "mean per-cell z on refinement", "1.066": "SD of that z",
    "12": "cells beyond two standard errors on refinement",
    "0.0009": "gap between the two type I error values",
    "504": "comparable cells in the refined grid",
    "16": "cells in the Delta_sig resolution block",
    "34": "bisection steps after refinement",
    "0.2": "the 0.2 pp movement threshold", "0.036": "mean movement of Delta_sig",
    "0.124": "largest movement of Delta_sig",
    "1.07": "largest movement as a percentage of Delta_tip",
    "1.6": "the 1.6 in 1.6 x 10^-9", "1.97": "the bracket reported after the fixed guard",
    "0.5": "0.5 per cent level of the one-sided criterion",
    "500": "500 patients per arm",
}


# ---------------------------------------------------------------------------
# Table 1
# ---------------------------------------------------------------------------
TABLE1_RE = re.compile(
    r"^\|\s*(0\.\d\d)\s*\|\s*(0\.\d\d)\s*\(p₁\s*=\s*(0\.\d+)\)\s*\|"
    r"\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|\s*$",
    re.M)


def check_table1(ms: str) -> None:
    rows = TABLE1_RE.findall(ms)
    guard("table1_shape", len(rows) == 9, f"{len(rows)} data rows, expected 9")
    if len(rows) != 9:
        return
    bad = []
    for p0s, ors, p1s, *cells in rows:
        p0, or_true = float(p0s), float(ors)
        p1 = legacy.p1_of(or_true, p0)
        if d(p1, 3) != p1s:
            bad.append(f"p1 label {p1s} != {d(p1, 3)} at p0={p0s} OR={ors}")
        for j, se in enumerate((0.70, 0.80, 0.90, 0.98)):
            want = d(legacy.delta_tip(or_true, p0, se, 0.0), 1)
            if cells[j] != want:
                bad.append(f"p0={p0s} OR={ors} se={se}: printed {cells[j]}, closed form {want}")
    guard("table1_every_cell_matches_the_closed_form", not bad,
          "all 36 cells and all 9 p1 labels reproduce" if not bad else "; ".join(bad[:4]))

    # the grid the table spans must be a subset of the certified evaluation grid
    spec = json.loads(SPEC.read_text(encoding="utf-8"))
    grid = spec["scenario_grid"]
    inside = (set(float(r[0]) for r in rows) <= set(grid["p0_true_reference_risk"])
              and set(float(r[1]) for r in rows) <= set(grid["or_true"])
              and {0.70, 0.80, 0.90, 0.98} <= set(grid["se_common"]))
    guard("table1_scenarios_lie_inside_the_certified_grid", inside,
          f"spec p0 {grid['p0_true_reference_risk']}, or {grid['or_true']}, se {grid['se_common']}")


# ---------------------------------------------------------------------------
# eMethods tables
# ---------------------------------------------------------------------------
def check_emethods(em: str, legacy_json: dict) -> None:
    p1_rows = legacy_json["part1_stated_model"]["rows"]
    bad = []
    for r in p1_rows:
        want = [d(r["delta_pp"], 2), d(r["analytic_observed_or"], 4),
                d(r["mc_median_observed_or"], 4), d(r["mc_ci2.5"], 3), d(r["mc_ci97.5"], 3)]
        star = r"\*{0,2}"
        pattern = re.compile(
            r"\|\s*" + star + re.escape(d(r["delta_multiple_of_tip"], 1)) + star + r"\s*\|"
            + r"\s*" + star + re.escape(want[0]) + star + r"\s*\|"
            + r"\s*" + star + re.escape(want[1]) + star + r"\s*\|"
            + r"\s*" + star + re.escape(want[2]) + star + r"\s*\|"
            + r"\s*" + star + re.escape(want[3]) + r"\s*[–-]\s*"
            + re.escape(want[4]) + star + r"\s*\|")
        if not pattern.search(em):
            bad.append(f"mult {r['delta_multiple_of_tip']}: expected {want}")
    guard("emethods_monte_carlo_table_matches_the_legacy_run", not bad,
          "all 5 rows x 5 values reproduce from reports/simulation_tipping_point.json"
          if not bad else "; ".join(bad[:3]))

    # part 2 table: effective (se, f) and observed risk per arm
    p2 = legacy_json["part2_boundary_model"]
    arms = [("Reference", p2["reference_arm"], legacy.GRADES_REF, 0.15, 0.10),
            ("Index (same process)", p2["index_arm_common_process"], legacy.GRADES_IDX, 0.15, 0.10)]
    bad = []
    for label, stored, dist, a, b in arms:
        se_eff, f_eff = legacy.effective(dist, a, b)
        risk = legacy.observed_risk_grade(dist, a, b)
        want = [d(a, 2), d(b, 2), d(se_eff, 4), d(f_eff, 4), d(risk, 4)]
        row = (rf"\|\s*{re.escape(label)}\s*\|\s*{want[0]}\s*\|\s*{want[1]}\s*\|"
               rf"\s*{want[2]}\s*\|\s*{want[3]}\s*\|\s*{want[4]}\s*\|")
        if not re.search(row, em):
            bad.append(f"{label}: printed row not found, expected {want}")
        if abs(se_eff - stored["se"]) > TOL or abs(f_eff - stored["f"]) > TOL \
                or abs(risk - stored["observed_risk"]) > TOL:
            bad.append(f"{label}: recomputed disagrees with the stored legacy record")
    guard("emethods_boundary_table_recomputes_from_the_model", not bad,
          "both arms reproduce from effective() and observed_risk_grade()" if not bad else "; ".join(bad[:3]))

    # OR = 1 contour table, all 9 rows x 3 printed values
    contour = legacy_json["part3_surface"]["or1_contour"]
    bad = []
    for c in contour:
        if c["a1_for_null"] is None:
            continue
        want = [d(c["b1"], 2), d(c["a1_for_null"], 4), d(c["delta_effective_pp"], 2),
                d(c["se_index"], 4)]
        bold = r"\*{0,2}"
        row = (rf"\|\s*{bold}{re.escape(want[0])}{bold}\s*\|\s*{bold}{re.escape(want[1])}{bold}\s*\|"
               rf"\s*{bold}{re.escape(want[2])}{bold}\s*\|\s*{bold}{re.escape(want[3])}{bold}\s*\|")
        if not re.search(row, em):
            bad.append(f"b1={c['b1']}: expected {want}")
    guard("emethods_or1_contour_table_matches_the_legacy_run", not bad,
          "all 9 contour rows reproduce" if not bad else "; ".join(bad[:3]))

    # grid size: the manuscript says 72 combinations; the grid is 16 x 9
    grid = legacy_json["part3_surface"]["grid"]
    n_a1 = len({r["a1"] for r in grid})
    n_b1 = len({r["b1"] for r in grid})
    claimed = re.search(r"steps of 0\.05 \((\d+) levels\) and b. from 0 to 0\.40 in steps of 0\.05 "
                        r"\((\d+) levels\), giving (\d+) combinations", em)
    if claimed:
        la, lb, tot = (int(x) for x in claimed.groups())
        detail = (f"printed {la} x {lb} = {tot}, actual {n_a1} x {n_b1} = {n_a1 * n_b1}")
        ok = (la, lb, tot) == (n_a1, n_b1, n_a1 * n_b1)
    else:
        single = re.search(r"\((\d+) combinations\)", em)
        ok = False
        detail = (f"the level counts are not stated; a bare combination count "
                  f"{single.group(1) if single else 'absent'} cannot be checked against the axes")
    guard("emethods_two_way_grid_size_is_correct", ok, detail)
    guard("emethods_one_way_slice_endpoints_match",
          d(legacy_json["part3_surface"]["one_way_slice"][0]["observed_or"], 4) == "0.5606"
          and d(legacy_json["part3_surface"]["one_way_slice"][-1]["observed_or"], 4) == "1.0813",
          "0.5606 at a1 = 0.05 and 1.0813 at a1 = 0.80")

    # shape robustness summary
    p4 = legacy_json["part4_shape_robustness"]
    lo, hi = p4["apparent_delta_abs_range_pp"]
    pct = max(abs(r["apparent_delta_as_pct_of_tip"]) for r in p4["rows"])
    guard("emethods_shape_artefact_summary_matches",
          f"| {d(lo,2)} pp |" in em.replace("**", "") or d(lo, 2) == "0.00",
          f"abs range {d(lo, 2)} to {d(hi, 2)} pp")
    guard("emethods_shape_artefact_maximum_is_four_five",
          d(hi, 2) == "4.50" and f"**{d(hi, 2)} pp**" in em,
          f"max {d(hi, 2)} pp")
    guard("emethods_shape_artefact_pct_of_threshold",
          d(pct, 1) == "23.2" and f"**{d(pct, 1)}%**" in em,
          f"max fraction {d(pct, 1)} per cent of Delta_tip")
    # the extreme really is at the two opposite corners of the concentration grid
    arg = sorted(((abs(r["apparent_delta_pp"]), r["conc_ref"], r["conc_idx"]) for r in p4["rows"]),
                 reverse=True)[:2]
    guard("emethods_shape_extremes_are_the_opposite_corners",
          {(a[1], a[2]) for a in arg} == {(0.30, 0.60), (0.60, 0.30)},
          f"extremes at {sorted((a[1], a[2]) for a in arg)}")
    zero_here = [r for r in p4["rows"] if r["conc_ref"] == r["conc_idx"]]
    guard("emethods_shape_artefact_vanishes_at_equal_concentration",
          len(zero_here) == 4 and all(abs(r["apparent_delta_pp"]) < 1e-12 for r in zero_here),
          f"{len(zero_here)} equal-concentration combinations, all zero")


# ---------------------------------------------------------------------------
# headline quantities
# ---------------------------------------------------------------------------
def check_headline(ms: str, em: str, base: dict, legacy_json: dict) -> None:
    pub = base["published_figures"]
    dt = legacy.delta_tip(0.60, 0.55, 0.90, 0.0)
    guard("headline_threshold_is_19_8",
          abs(dt - 19.8) < 0.005 and ms.count("19.8") >= 6,
          f"delta_tip(0.60, 0.55, 0.90, 0) = {dt:.4f} pp, printed {ms.count('19.8')} times")

    scan = [(se, legacy.delta_tip(0.60, 0.55, se, 0.0)) for se in (0.70, 0.98)]
    guard("se_scan_range_is_15_4_to_21_6",
          d(min(v for _, v in scan), 1) == "15.4" and d(max(v for _, v in scan), 1) == "21.6"
          and "15.4–21.6" in ms,
          ", ".join(f"se={se}: {v:.4f} -> printed {d(v, 1)}" for se, v in scan))

    contour = [c["delta_effective_pp"] for c in legacy_json["part3_surface"]["or1_contour"]
               if c["a1_for_null"] is not None]
    guard("two_way_family_range_is_17_70_to_26_02",
          d(min(contour), 2) == "17.70" and d(max(contour), 2) == "26.02",
          f"{min(contour):.4f} to {max(contour):.4f} pp, printed {d(min(contour), 2)} to {d(max(contour), 2)}")
    guard("manuscript_two_way_range_is_a_correct_rounding",
          "17.7 to 26.0 percentage points" in ms,
          "26.02 printed to one decimal as 26.0; 17.70 to 17.7")

    a1 = legacy_json["part2_boundary_model"]["index_upgrade_rate_for_null"]
    guard("index_upgrade_rate_for_a_null_is_71_percent",
          abs(a1 - 0.7119815668202765) < TOL and a1 * 100 > 70.5 and "71%" in ms,
          f"a1 = {a1:.4f}")
    guard("reference_upgrade_rate_is_15_percent",
          abs(legacy_json["part2_boundary_model"]["reference_arm"]["a"] - 0.15) < TOL
          and "against 15% in the reference arm" in ms,
          "a0 = 0.15")

    p2 = legacy_json["part2_boundary_model"]
    ratio = p2["ratio_needed_to_closed_form"]
    guard("closed_form_to_boundary_ratio_is_0_997",
          abs(ratio - 0.997) < 0.0005 and abs(p2["delta_tip_closed_form_pp"] - 19.83) < 0.005
          and abs(p2["delta_needed_for_null_pp"] - 19.78) < 0.005,
          f"{p2['delta_tip_closed_form_pp']:.4f} vs {p2['delta_needed_for_null_pp']:.4f}, ratio {ratio:.6f}")
    guard("index_effective_false_positive_at_crossing_is_0_2678",
          abs(p2["index_effective_f_at_null"] - 0.2678) < 0.00005,
          f"{p2['index_effective_f_at_null']:.6f}")

    # observed-for-true substitution
    p0, p1 = legacy.P0, legacy.p1_of(0.60, 0.55)
    q0 = legacy.observed_risk_effective(p0, 0.90, 0.0)
    q1 = legacy.observed_risk_effective(p1, 0.90, 0.0)
    naive = legacy.delta_tip_from_risks(0.90, 0.0, q0, q1) if hasattr(legacy, "delta_tip_from_risks") \
        else 100.0 * (q0 - q1) * 0.90 / (1 - q1)
    under = 100.0 * (dt - naive) / dt
    guard("observed_substitution_is_16_6_pp_and_16_percent",
          abs(naive - 16.6) < 0.05 and abs(under - 16.0) < 0.6
          and "16.6 percentage points" in ms and "underestimate of 16%" in ms,
          f"observed risks {q0:.4f}/{q1:.4f} give {naive:.2f} pp, {under:.2f} per cent low")
    guard("observed_risks_are_0_495_and_0_381",
          d(q0, 3) == "0.495" and d(q1, 3) == "0.381"
          and "the observed risks are 0.495 and 0.381" in ms,
          f"{q0:.6f}, {q1:.6f}")

    # same-denominator conversions of the two external disagreement figures
    conv = {}
    for share in (0.106, 0.258):
        vals = []
        for base_risk in (0.50, 0.55, 0.60):
            p1x = legacy.p1_of(0.60, base_risk)
            vals.append(100.0 * (share / 2) / (1 - p1x))
        conv[share] = vals
    guard("same_denominator_conversion_of_10_6_percent",
          abs(min(conv[0.106]) - 8.5) < 0.05 and abs(max(conv[0.106]) - 10.1) < 0.05
          and "8.5 to 10.1 percentage points" in ms,
          f"{min(conv[0.106]):.2f} to {max(conv[0.106]):.2f} pp")
    guard("same_denominator_conversion_of_25_8_percent",
          abs(min(conv[0.258]) - 20.6) < 0.05 and abs(max(conv[0.258]) - 24.5) < 0.05
          and "20.6 to 24.5" in ms,
          f"{min(conv[0.258]):.2f} to {max(conv[0.258]):.2f} pp")

    # f = 0 is the maximum of the family
    fam = [(legacy.delta_tip(0.60, 0.55, 0.90, f), f) for f in (0.0, 0.05, 0.10, 0.15)]
    guard("f_zero_is_the_largest_in_the_family",
          fam[0][0] == max(v for v, _ in fam) and fam[0][1] == 0.0
          and "the value at f = 0 is the largest in this family" in ms,
          ", ".join(f"f={f}: {v:.2f}" for v, f in fam))

    # non-differential error preserves the sign and attenuates
    ratios = []
    for or_true in (0.60, 0.70, 0.80):
        p1x = legacy.p1_of(or_true, 0.55)
        r0 = legacy.observed_risk_effective(0.55, 0.90, 0.05)
        r1 = legacy.observed_risk_effective(p1x, 0.90, 0.05)
        ratios.append(legacy.or_from_risks(r0, r1))
    guard("non_differential_error_attenuates_but_never_reverses",
          all(0 < r < 1 for r in ratios) and all(abs(0.55 - 0.423) * (0.90 - 0.05) > 0 for _ in [0]),
          ", ".join(f"OR_obs {r:.4f}" for r in ratios))

    # the baseline step must have re-derived all of these already
    guard("baseline_artefact_agrees_with_this_run",
          pub["main_delta_tip_pp"] == 19.8 and pub["ratio_to_closed_form"] == 0.997
          and pub["common_process_artefact_pp"] == 4.5,
          "baseline published_figures matches the live recomputation")
    guard("baseline_artefact_declares_what_it_did_not_recompute",
          len(base["not_recomputed_here"]) >= 5,
          f"{len(base['not_recomputed_here'])} external inputs listed as not recomputed")


# ---------------------------------------------------------------------------
# C. claims v12 must not make, and findings it does not yet contain
# ---------------------------------------------------------------------------
# Pattern guards for claims that were once made and were removed. Each must be
# non-vacuous, so each has a matching injection in the battery below.
FORBIDDEN_MS = [
    ("ci_threshold", re.compile(r"(?:19\.8|Δ_tip|threshold)[^.]{0,80}confidence interval[^.]{0,60}(?:becomes|null|cross)"),
     "Delta_tip is the point at which the point estimate is null"),
    ("significance_threshold", re.compile(r"(?:19\.8|Δ_tip)[^.]{0,50}\bis\b[^.]{0,30}significan"),
     "Delta_tip carries no significance claim"),
    ("one_in_five", re.compile(r"one in five"),
     "the share-of-assignments misstatement was removed"),
    ("observed_as_097", re.compile(r"observed as 0\.97"),
     "the OR = 0.97 misstatement was removed"),
]

# "significan*" is legal in exactly one sentence of v12 -- the caveat that the
# single published between-arm comparison returned a non-significant result --
# and illegal everywhere else. Whitelisting the phrase is robust; a proximity
# window around "threshold" is not, because the distance depends on how the
# sentence is punctuated.
SIGNIFICANCE_ALLOWED = (
    "non-significant result does not establish equivalence",   # the one published comparison
    "loses significance",                                     # Delta_sig, twinned with Delta_tip
    "criterion for significance",                             # the warning against using Delta_tip as one
)

# "smallest" / "minimum" / "lower bound" are legitimate words here, but only in
# negated form, and the negation has to sit in the same CLAUSE. A sentence-level
# check is not enough: the manuscript's own disclaimer ("the one-way value is
# not the smallest in it") would excuse a violation asserted in the same
# sentence, which is exactly what the mutation below does.
NEGATION = re.compile(r"\b(?:not|no|never|cannot|does not|do not|declin\w*)\b", re.I)
MINIMALITY = re.compile(r"\b(?:smallest|minimum|lower bound)\b", re.I)
CLAUSE_SPLIT = re.compile(r"[;:.!?]|,\s+(?:and|but|so|which|while|where)\b|\s+and\s+")
# "minimum important difference" is a named clinical concept, not a claim about
# the threshold, so it is removed before the minimality scan.
MID = re.compile(r"minimum (?:clinically )?important difference", re.I)
# "minimum strength" and "minimum-strength-to-explain-away" are the name of an
# existing parameter (VanderWeele & Li 2019), not a claim about this threshold.
NAMED_CONSTRUCT = re.compile(
    r"minimum-strength-to-explain-away|minimum strength(?: of differential measurement error)?",
    re.I)


def check_overclaims(ms: str) -> None:
    hits = []
    for name, pattern, why in FORBIDDEN_MS:
        m = pattern.search(ms)
        if m:
            hits.append(f"{name}: {m.group(0)[:70]!r} ({why})")

    # every occurrence of "significan*" must be the whitelisted caveat
    scrubbed = ms
    for phrase in SIGNIFICANCE_ALLOWED:
        scrubbed = scrubbed.replace(phrase, "")
    stray = re.findall(r"\w*significan\w*", scrubbed, re.I)
    if stray:
        hits.append(f"significance language outside the one permitted caveat: {stray[:4]}")
    guard("manuscript_makes_no_claim_stronger_than_it_can_support", not hits,
          f"{len(FORBIDDEN_MS)} forbidden patterns absent; 'significan*' occurs only in "
          f"the permitted caveat" if not hits else "; ".join(hits))

    scan = NAMED_CONSTRUCT.sub("named-construct", MID.sub("mid", ms))
    unnegated = []
    for clause in CLAUSE_SPLIT.split(scan):
        if MINIMALITY.search(clause) and not NEGATION.search(clause):
            unnegated.append(clause.strip()[:90])
    guard("every_minimality_word_is_negated_in_its_own_clause", not unnegated,
          "each occurrence of smallest/minimum/lower bound sits in a negated clause"
          if not unnegated else "; ".join(unnegated[:3]))

    # conflating the apparent asymmetry with the closed-form failure is the
    # specific error step 6 was written to prevent
    conflating = re.search(
        r"closed form is (?:wrong|off|inaccurate) by (?:up to )?[\d.]+ percentage points",
        ms)
    guard("manuscript_does_not_conflate_apparent_asymmetry_with_closed_form_failure",
          conflating is None,
          "no such sentence" if conflating is None else conflating.group(0))


def check_w1_present(cal: dict, dsig: dict, cont: dict, ms: str, em: str) -> None:
    pend = []

    pend.append({
        "finding": "the closed form fails by -2.13 to +6.19 pp across the layer B grid",
        "value": [round(min(r["closed_form_failure_pp"] for r in cont["rows"]), 2),
                  round(max(r["closed_form_failure_pp"] for r in cont["rows"]), 2)],
        "rows": len(cont["rows"]),
        "manuscript_status": "present in Results, How the threshold behaves under sampling",
        "source": "reports/estimator_eval_null_contour.json",
    })
    pend.append({
        "finding": "5 of 73 layer B combinations have no reachable OR = 1 contour",
        "value": cont["combinations"]["without_a_contour"],
        "rows": cont["combinations"]["total"],
        "manuscript_status": "present in Results",
        "source": "reports/estimator_eval_null_contour.json",
    })
    pend.append({
        "finding": "pooled type I error at the true null",
        "value": round(cal["summary"]["pooled_type_I_error"], 4),
        "rows": cal["summary"]["type_I_cells"],
        "manuscript_status": "present in Results and the abstract",
        "source": "reports/estimator_eval_calibration.json",
    })
    pend.append({
        "finding": "mean Delta_sig / Delta_tip",
        "value": round(dsig["summary"]["delta_sig_over_delta_tip_mean"], 4),
        "rows": dsig["summary"]["delta_sig_cells"],
        "manuscript_status": "present in Results",
        "source": "reports/estimator_eval_delta_sig_coverage.json",
    })
    pend.append({
        "finding": "two-layer coverage, uncorrected then corrected",
        "value": [round(dsig["block2_coverage"]["pooled_coverage_layer1_uncorrected"], 4),
                  round(dsig["block2_coverage"]["pooled_coverage_layer2_corrected"], 4)],
        "rows": dsig["summary"]["coverage_cells"],
        "manuscript_status": "present in Results",
        "source": "reports/estimator_eval_delta_sig_coverage.json",
    })

    # Every one of the five W1 findings must now be *present* in the submitted
    # text, each in the file where it belongs. A finding that is absent is a
    # failure, not a note: bmcmrm_v1 claims these results in its abstract.
    lo, hi = pend[0]["value"]
    absent = []
    checks = [
        ("layer-B failure range in the manuscript",
         f"−{abs(lo):.2f} to +{hi:.2f}" in ms, ms),
        ("no-contour count in the manuscript", "the remaining 5" in ms, ms),
        ("pooled type I error in the manuscript", "0.0525" in ms, ms),
        ("Delta_sig shortfall in the manuscript", "1.91 percentage points" in ms, ms),
        ("Delta_sig/Delta_tip mean in the eMethods", "0.8583" in em, em),
        ("uncorrected coverage in the manuscript", "87.8%" in ms, ms),
        ("corrected coverage in the manuscript", "95.2%" in ms, ms),
    ]
    for label, ok, _blob in checks:
        if not ok:
            absent.append(label)
    guard("w1_findings_are_present_where_they_are_claimed", not absent,
          f"{len(checks)} presence checks pass" if not absent else "; ".join(absent))

    note("W1 findings and where bmcmrm_v1 carries each (presence asserted above):")
    for p in pend:
        note(f"  - {p['finding']} = {p['value']}  [{p['rows']} cells]  ({p['manuscript_status']})")
    _lines.append(json.dumps({"w1_findings": pend}, ensure_ascii=False, indent=1))

    statuses = {}
    for path in (BASE, DGP, CAL, DSIG):
        blob = json.loads(path.read_text(encoding="utf-8"))
        statuses[path.stem.replace("estimator_eval_", "")] = bool(blob["summary"]["all_passed"])
    statuses["null_contour"] = bool(cont["checks_passed"])
    design = json.loads(SPEC.read_text(encoding="utf-8"))
    statuses["design"] = all(c["passed"] for c in design["self_checks"])
    guard("all_six_step_artefacts_report_success", all(statuses.values()),
          ", ".join(f"{k}={v}" for k, v in statuses.items()))


# ---------------------------------------------------------------------------
# D. sampling behaviour, recomputed from the artefacts
# ---------------------------------------------------------------------------
def check_sampling_recompute(ms: str, em: str, cal: dict, dsig: dict,
                             cont: dict, conv: dict) -> None:
    """Recompute every sampling-behaviour number from the JSON that produced it.

    The accounting guard proves a number was *declared*; this guard proves it is
    the number the artefact holds. Both are needed: a declared value that has
    drifted from its source passes accounting and fails here.
    """
    checks: list[tuple[str, bool, str]] = []
    pooled = cal["tables"]["I_type_I_error"]["pooled"]
    rate, mcse = pooled["rate"], pooled["mc"]
    checks += [
        # the abstract AND the Results both state this; a mutation that strips
        # one occurrence must fail, so the presence check is a count, not a
        # substring test (this is the W7 lesson: same phrase, two sites).
        ("pooled type I error, manuscript",
         d(rate, 4) == "0.0525" and ms.count("0.0525") >= 2,
         f"{rate!r} -> {d(rate, 4)}, {ms.count('0.0525')} occurrences"),
        ("pooled type I error, eMethods", d(rate, 5) == "0.05252" and "0.05252" in em,
         f"{d(rate, 5)}"),
        ("Monte Carlo SE, manuscript",
         d(mcse, 4) == "0.0008" and ms.count("0.0008") >= 2,
         f"{mcse!r} -> {d(mcse, 4)}, {ms.count('0.0008')} occurrences"),
        ("Monte Carlo SE, eMethods", d(mcse, 5) == "0.00077" and "0.00077" in em,
         f"{d(mcse, 5)}"),
        ("type I error cell count", str(cal["summary"]["type_I_cells"]) in ms,
         f"{cal['summary']['type_I_cells']} cells"),
    ]
    rates = [c["reject_rate"] for c in cal["tables"]["I_type_I_error"]["cells"]]
    checks += [
        ("single-cell type I error range",
         d(min(rates), 4) == "0.0370" and d(max(rates), 4) == "0.0665"
         and "0.0370" in em and "0.0665" in em,
         f"{d(min(rates), 4)} to {d(max(rates), 4)}"),
    ]

    b1 = dsig["block1_delta_sig"]
    ratio, cells = b1["ratio_overall"]["mean"], dsig["summary"]["delta_sig_cells"]
    short = [c["delta_sig_pp"] - c["delta_tip_pp"] for c in b1["cells"]]
    checks += [
        ("Delta_sig cell count", str(cells) in ms and str(cells) in em, f"{cells} cells"),
        ("mean Delta_sig/Delta_tip", d(ratio, 4) == "0.8583" and "0.8583" in em,
         f"{ratio!r} -> {d(ratio, 4)}"),
        ("mean shortfall, manuscript", d(abs(sum(short) / len(short)), 2) == "1.91"
         and ms.count("1.91") >= 2 and "1.91 percentage points" in ms,
         f"{sum(short) / len(short):.4f} pp, {ms.count('1.91')} occurrences"),
        ("largest shortfall, eMethods", d(abs(min(short)), 2) == "5.63" and "5.63" in em,
         f"{min(short):.4f} pp"),
        ("no cell reaches Delta_tip", all(s < 0 for s in short),
         f"{sum(1 for s in short if s >= 0)} of {len(short)} cells at or above Delta_tip"),
    ]
    for n, want in (("100", "0.7219"), ("500", "0.8700"), ("1000", "0.9070"), ("2000", "0.9341")):
        got = d(b1["ratio_by_n"][n]["mean"], 4)
        checks.append((f"mean ratio at n={n}", got == want and want in em, f"{got}"))

    cov = dsig["block2_coverage"]
    l1, l2 = cov["pooled_coverage_layer1_uncorrected"], cov["pooled_coverage_layer2_corrected"]
    checks += [
        ("coverage cell count", str(dsig["summary"]["coverage_cells"]) in ms
         and str(dsig["summary"]["coverage_cells"]) in em, f"{dsig['summary']['coverage_cells']}"),
        ("uncorrected coverage", d(l1, 4) == "0.8779" and "0.8779" in em
         and d(l1 * 100, 1) == "87.8" and "87.8%" in ms, f"{l1!r}"),
        ("corrected coverage", d(l2, 4) == "0.9516" and "0.9516" in em
         and d(l2 * 100, 1) == "95.2" and "95.2%" in ms, f"{l2!r}"),
    ]
    for n, want_pct, want_em in (("100", "94.2", "0.9420"), ("2000", "79.6", "0.7959")):
        rows = [c for c in cov["cells"] if str(c["n_per_arm"]) == n]
        got = sum(c["coverage_layer1_uncorrected"] for c in rows) / len(rows)
        checks.append((f"uncorrected coverage at n={n}",
                       d(got, 4) == want_em and want_em in em
                       and d(got * 100, 1) == want_pct and f"{want_pct}%" in ms,
                       f"{got!r} -> {d(got, 4)} over {len(rows)} cells"))

    fails = [r["closed_form_failure_pp"] for r in cont["rows"]]
    checks += [
        ("contour combination counts",
         str(cont["combinations"]["with_a_contour"]) in ms
         and str(cont["combinations"]["total"]) in ms,
         f"{cont['combinations']}"),
        ("layer-B failure range, manuscript",
         d(abs(min(fails)), 2) == "2.13"
         and ms.count(f"−{d(abs(min(fails)), 2)} to +{d(max(fails), 2)}") >= 2,
         f"{min(fails):.4f} to {max(fails):.4f}"),
        ("layer-B failure range, eMethods",
         d(max(fails), 2) == "6.19" and "6.19" in em, f"{d(max(fails), 2)}"),
        ("failures vanish exactly at equal sensitivity",
         sum(1 for v in fails if abs(v) < 1e-9) == 16, 
         f"{sum(1 for v in fails if abs(v) < 1e-9)} of {len(fails)} are zero"),
    ]

    h = conv["summary"]["headline"]
    checks += [
        ("paired z on doubling, mean",
         d(h["mean_z"], 4) == "0.0230" and "0.0230" in em, f"{h['mean_z']!r}"),
        ("paired z on doubling, SD",
         d(h["sd_z"], 4) == "1.0046" and "1.0046" in em, f"{h['sd_z']!r}"),
        ("pooled type I, coarse grid",
         d(h["pooled_type_I_coarse"], 6) == "0.051678" and "0.051678" in em,
         f"{h['pooled_type_I_coarse']!r}"),
        ("pooled type I, refined grid",
         d(h["pooled_type_I_fine"], 6) == "0.051600" and "0.051600" in em,
         f"{h['pooled_type_I_fine']!r}"),
        ("grid refinement difference",
         d(abs(h["pooled_type_I_difference"]), 6) == "0.000078" and "0.000078" in em,
         f"{h['pooled_type_I_difference']!r}"),
        ("largest Delta_sig movement",
         d(h["max_delta_sig_shift_pp"], 3) == "0.124" and "0.124" in em,
         f"{h['max_delta_sig_shift_pp']!r}"),
        ("grid sizes", str(h["coarse_cells"]) in ms and str(h["fine_cells"]) in ms
         and f"{h['coarse_cells']} to {h['fine_cells']} cells" in ms,
         f"{h['coarse_cells']} -> {h['fine_cells']}"),
    ]

    bad = [f"{n}: {detail}" for n, ok, detail in checks if not ok]
    guard("sampling_behaviour_numbers_recompute_from_their_artefacts", not bad,
          f"{len(checks)} values recomputed, all match" if not bad else "; ".join(bad[:5]))


# ---------------------------------------------------------------------------
# eTable 1
# ---------------------------------------------------------------------------
def check_etable1(ms: str, et: str) -> None:
    rows = re.findall(r"^\|\s*(S\d)\s*\|[^|]*\|[^|]*\|\s*([\d,]+)\s*\|", et, re.M)
    s1 = {k: v for k, v in rows}
    guard("etable1_s1_yield_matches_the_manuscript",
          s1.get("S1") == "1,188" and "1,188" in ms,
          f"eTable 1 S1 = {s1.get('S1')}, manuscript says 1,188")
    combined = re.search(r"union of S1–S8 \|\s*\*\*([\d,]+)\*\*", et)
    guard("etable1_records_all_its_own_yields",
          combined is not None and len(rows) == 8,
          f"{len(rows)} query rows, combined total {combined.group(1) if combined else 'absent'}")
    g = re.findall(r"^\|\s*(G\d)\s*\|[^|]*\|[^|]*\|\s*\*\*(\d+)\*\*", et, re.M)
    guard("etable1_gap_searches_are_recorded_with_their_yields",
          len(g) == 3 and dict(g) == {"G1": "1", "G2": "0", "G3": "0"},
          f"{dict(g)}")


# ---------------------------------------------------------------------------
# abstract: the one hard length gate, and the two copies of it
# ---------------------------------------------------------------------------
ABSTRACT_LABELS = ("Background.", "Methods.", "Results.", "Conclusions.")


def _abstract_body(ms: str, ab: str) -> tuple[str, str]:
    """The structured abstract, normalised, as carried by each document.

    Three traps this exists to close.

    1. The standalone copy in abstract.md is prefixed with title / article-type
       / journal lines, so a whitespace count of the whole file reads ~414 words
       and would look like a violation.
    2. The manuscript body is delimited by the next Markdown heading, which is
       why the section has to be cut rather than counted in place.
    3. The two copies carry the section labels differently -- the manuscript
       keeps "Background." inline in the section body while the standalone copy
       opens with it -- so the standalone body is re-prefixed here.  Without
       that, the comparison reports a one-word drift that is purely an artefact
       of the markup, and a guard that cries wolf is a guard that gets ignored.
    """
    norm = lambda s: re.sub(r"\s+", " ", s.replace("**", "")).strip()  # noqa: E731
    m = re.search(r"^##\s*Abstract\s*$(.+?)(?=^##\s)", ms, re.S | re.M)
    ms_raw = norm(m.group(1)) if m else ""
    m = re.search(r"\*\*Background\.\*\*(.+?)\*\*Keywords:\*\*", ab, re.S)
    ab_raw = norm("Background. " + m.group(1)) if m else ""
    return ms_raw, ab_raw


def _abstract_prose(body: str) -> str:
    """The abstract with its four structural labels removed.

    This is what has to match word for word; the labels are presentation.
    """
    for lab in ABSTRACT_LABELS:
        body = re.sub(r"(?<!\S)" + re.escape(lab) + r"(?!\S)", "", body)
    return re.sub(r"\s+", " ", body).strip()


def check_abstract_limit(ms: str, ab: str, cp: str) -> None:
    ms_raw, ab_raw = _abstract_body(ms, ab)
    n_ms, n_ab = len(ms_raw.split()), len(ab_raw.split())
    prose = _abstract_prose(ms_raw)

    # the limit is asserted on the label-inclusive count, i.e. the worst case
    guard("abstract_stays_within_the_journal_limit",
          bool(ms_raw) and n_ms <= ABSTRACT_LIMIT and n_ab <= ABSTRACT_LIMIT,
          f"manuscript {n_ms} words incl. labels ({len(prose.split())} of prose), "
          f"standalone {n_ab} words, limit {ABSTRACT_LIMIT}")

    # the submitted copy and the standalone copy must be the same text; they are
    # entered into two different boxes of the submission system
    ab_prose = _abstract_prose(ab_raw)
    guard("abstract_copy_is_identical_to_the_manuscript",
          bool(ab_prose) and prose == ab_prose,
          "standalone copy matches the manuscript verbatim"
          if ab_prose and prose == ab_prose
          else f"the two copies differ ({len(prose.split())} vs {len(ab_prose.split())} words of prose)")

    # A single fact stated in two places drifts.  Both declarations are read and
    # required to agree with each other *and* with the measured count, so that
    # editing one of the two is caught rather than merely made inconsistent.
    # Only lines about the abstract are scanned: a bare "300 words" also occurs
    # in the figure-legend row and must not be mistaken for a declaration here.
    decls = [n for ln in cp.splitlines() if "abstract" in ln.lower()
             for n in re.findall(r"(\d{3})\s*(?:words?\b|at last check)", ln, re.I)]
    guard("compliance_states_the_same_abstract_word_count",
          bool(decls) and set(decls) == {str(n_ms)},
          f"reporting_compliance.md declares {sorted(set(decls)) or 'nothing'}, "
          f"measured {n_ms}")


# ---------------------------------------------------------------------------
# accounting
# ---------------------------------------------------------------------------
def check_accounting(text: str, tables: list[tuple[dict, str]]) -> None:
    seen = tokens(text)
    unknown = []
    for tok, lines in sorted(seen.items(), key=lambda kv: (len(kv[0]), kv[0])):
        if any(tok in t for t, _ in tables):
            continue
        unknown.append(f"{tok} (lines {sorted(lines)})")
    known = sum(1 for tok in seen if any(tok in t for t, _ in tables))
    guard("every_numeric_token_is_accounted_for", not unknown,
          f"{len(seen)} distinct tokens, all classified" if not unknown
          else f"{len(unknown)} unclassified: " + "; ".join(unknown[:8]))


# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------
def run() -> bool:
    texts = load_texts()
    legacy_json = json.loads(LEGACY_JSON.read_text(encoding="utf-8"))
    base = json.loads(BASE.read_text(encoding="utf-8"))
    cal = json.loads(CAL.read_text(encoding="utf-8"))
    dsig = json.loads(DSIG.read_text(encoding="utf-8"))
    cont = json.loads(CONT.read_text(encoding="utf-8"))
    conv = json.loads(CONV.read_text(encoding="utf-8"))

    ms_body = texts["ms"][: texts["ms"].index("## References")]

    check_table1(texts["ms"])
    check_emethods(texts["em"], legacy_json)
    check_headline(texts["ms"], texts["em"], base, legacy_json)
    check_overclaims(ms_body)
    check_w1_present(cal, dsig, cont, texts["ms"], texts["em"])
    check_sampling_recompute(texts["ms"], texts["em"], cal, dsig, cont, conv)
    check_etable1(texts["ms"], texts["et"])
    check_abstract_limit(texts["ms"], texts["ab"], texts["cp"])
    check_accounting(ms_body, [(STRUCTURAL_MS, "structural"), (EXTERNAL_MS, "external"),
                               (PARAM_MS, "parameter"), (DERIVED_MS, "derived")])
    check_accounting(texts["em"], [(STRUCTURAL_EM, "structural"), (EXTERNAL_EM, "external"),
                                   (PARAM_EM, "parameter"), (DERIVED_EM, "derived")])

    failed = [n for n, ok, _ in _guards if not ok]
    _lines.append("")
    _lines.append(f"guards: {len(_guards)}, failed: {len(failed)}")
    for n in failed:
        _lines.append(f"  FAILED: {n}")
    REPORT.write_text("\n".join(_lines), encoding="utf-8")
    print(f"\nreport -> {REPORT}")
    return not failed


def selftest() -> bool:
    """Every injection must produce a failing guard that the clean run does not.

    Two failure modes are separated deliberately, because both have happened in
    this project before:

      * an anchor that is not present makes the loader raise SystemExit, the
        process dies with a non-zero code, and a naive "did it fail?" test calls
        that a successful catch. The guard here is on the *set* of failing
        guards, not on the exit code, and the loader is required to have
        rewritten the file (checked by asserting the anchor was consumed).
      * the clean run currently has one real failure, so "the run failed" is not
        evidence of anything. Only a guard absent from the baseline set counts.
    """
    import subprocess

    script = Path(__file__).resolve()

    # every subprocess is routed away from the delivery tree, so the shipped
    # report always describes the clean run and never an injected one
    import tempfile
    scratch = tempfile.mkdtemp(prefix="vdr_selftest_")

    # what the delivery report said before the battery; the battery must not
    # change it, or the shipped summary would describe an injected run
    shipped_before = REPORT.read_text(encoding="utf-8") if REPORT.exists() else None

    def failing_guards(env_extra: dict[str, str]) -> tuple[set[str], str]:
        env = dict(os.environ, VDR_OUTDIR=scratch, **env_extra)
        p = subprocess.run([sys.executable, str(script)], env=env,
                           capture_output=True, text=True)
        names = {ln.split("]")[1].split("->")[0].strip()
                 for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")}
        return names, p.stderr

    baseline, _ = failing_guards({})
    print(f"baseline failing guards ({len(baseline)}): {sorted(baseline)}\n")
    print("injection battery: each mutation must add a NEW failing guard\n")

    survivors, invalid = [], []
    for name, (key, old, new, desc) in INJECTIONS.items():
        # the anchor must exist verbatim, or the injection proves nothing
        probe = {"ms": MS, "em": EM, "et": ET, "ab": AB, "cp": CP}[key].read_text(encoding="utf-8")
        if old not in probe:
            print(f"  INVALID   {name:22s} anchor absent from {key}")
            invalid.append(name)
            continue
        if old == new:
            print(f"  INVALID   {name:22s} mutation is a no-op")
            invalid.append(name)
            continue

        after, err = failing_guards({"VMN_INJECT": name})
        if "INJECTION ERROR" in err:
            print(f"  INVALID   {name:22s} {err.strip().splitlines()[-1][:70]}")
            invalid.append(name)
            continue
        new_names = after - baseline
        if not new_names:
            print(f"  SURVIVED  {name:22s} {desc}")
            survivors.append(name)
            continue
        print(f"  CAUGHT    {name:22s} {desc}")
        for n in sorted(new_names):
            print(f"              -> {n}")

    print()
    ok = True
    if invalid:
        print(f"INVALID: {len(invalid)} injections did not actually mutate anything: {invalid}")
        ok = False
    if survivors:
        print(f"INERT: {len(survivors)} injections added no new failing guard: {survivors}")
        ok = False

    shipped_after = REPORT.read_text(encoding="utf-8") if REPORT.exists() else None
    if shipped_after != shipped_before:
        print("LEAK: the injection battery wrote into the delivery report; "
              "the shipped summary now describes an injected run")
        ok = False

    import shutil
    shutil.rmtree(scratch, ignore_errors=True)

    if ok:
        print(f"all {len(INJECTIONS)} injections applied and caught")
    return ok


def main() -> int:
    if "--list" in sys.argv:
        for k, (_f, _o, _n, desc) in INJECTIONS.items():
            print(f"{k:24s} {desc}")
        return 0
    if "--selftest" in sys.argv:
        return 0 if selftest() else 1
    return 0 if run() else 1


if __name__ == "__main__":
    raise SystemExit(main())