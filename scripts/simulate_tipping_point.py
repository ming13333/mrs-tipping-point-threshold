# -*- coding: utf-8 -*-
"""Validation of the one-way tipping point for differential misclassification.

Three questions, in increasing order of how much they threaten the manuscript.

  PART 1 -- Is the closed form right?
      Monte-Carlo under the model the manuscript states (equal effective
      sensitivity, false-positive rate f in the reference arm and f + Delta in
      the index arm). Check that the empirical observed OR crosses 1 at
      Delta_tip. This is a numerical check of the algebra, not a new result.

  PART 2 -- Does it survive a more realistic error structure?
      The manuscript's model is stated in terms of an *effective* sensitivity
      and false-positive rate. Real disagreement concentrates at the mRS 2/3
      boundary. So model the error where it actually occurs: a true mRS 2 graded
      3 (upgrade, rate a) and a true mRS 3 graded 2 (downgrade, rate b), applied
      to arm-specific grade distributions. The effective (se, f) then follow
      from (a, b) and the distributions. Question: how does the asymmetry
      needed for a null, expressed at the level where error occurs, relate to
      the Delta_tip the closed form reports?

  PART 3 -- Can a NON-differential process look differential?
      This is the sharpest question. If both arms share the same (a, b) but have
      different grade distributions -- which they do, because the outcome
      differs -- then the effective se and f differ between arms even though the
      measurement process is identical. Quantify the size and direction of the
      resulting artefact, and check whether it can move the observed OR.

Run:  python scripts/simulate_tipping_point.py
Out:  reports/simulation_tipping_point.json
"""
import io
import json
import os
import random

random.seed(20260926)

# ---------------------------------------------------------------------------
# Grade distributions (illustrative, and labelled as such in the manuscript)
# ---------------------------------------------------------------------------
# Chosen so that the reference arm has a 55% poor-outcome risk and the index arm
# 42.3%, i.e. a true odds ratio of 0.60 at p0 = 0.55 -- the scenario tabulated in
# the manuscript. The shape within the favourable grades deliberately differs
# between arms, because a treatment that shifts patients out of the poor
# categories does not shift them proportionally within the good ones.
GRADES_REF = {0: 0.09, 1: 0.15, 2: 0.21, 3: 0.16, 4: 0.15, 5: 0.11, 6: 0.13}
GRADES_IDX = {0: 0.16, 1: 0.20, 2: 0.217, 3: 0.12, 4: 0.12, 5: 0.09, 6: 0.093}

P0 = sum(v for g, v in GRADES_REF.items() if g >= 3)   # 0.55
P1 = sum(v for g, v in GRADES_IDX.items() if g >= 3)   # 0.423

assert abs(P0 - 0.55) < 1e-9, P0
assert abs(P1 - 0.423) < 1e-9, P1
TRUE_OR = (P1 / (1 - P1)) / (P0 / (1 - P0))


def effective(dist, a, b):
    """Effective (sensitivity, false-positive rate) for the mRS 3-6 dichotomy.

    a = P(graded 3 | true 2)   -- upgrade, drives false positives
    b = P(graded 2 | true 3)   -- downgrade, drives false negatives
    Deaths and grades other than 2/3 are assumed not to move.
    """
    poor = sum(v for g, v in dist.items() if g >= 3)
    fav = 1 - poor
    p2 = dist[2]
    p3 = dist[3]
    f = a * p2 / fav if fav > 0 else 0.0          # P(observed >=3 | true <=2)
    se = 1 - b * p3 / poor if poor > 0 else 1.0   # P(observed >=3 | true >=3)
    return se, f


def observed_risk_grade(dist, a, b):
    """Exact P(observed >=3) under the boundary-error model."""
    poor = sum(v for g, v in dist.items() if g >= 3)
    return poor - b * dist[3] + a * dist[2]


def observed_risk_effective(p, se, f):
    """The manuscript's model, for cross-checking."""
    return p * se + (1 - p) * f


def odds(risk):
    return risk / (1 - risk)


def or_from_risks(r0, r1):
    return odds(r1) / odds(r0)


# ---------------------------------------------------------------------------
# The closed form from the manuscript
# ---------------------------------------------------------------------------
def p1_of(or_true, p0):
    o = or_true * p0 / (1 - p0)
    return o / (1 + o)


def delta_tip(or_true, p0, se, f=0.0):
    p1 = p1_of(or_true, p0)
    return 100.0 * (p0 - p1) * (se - f) / (1 - p1)


# ---------------------------------------------------------------------------
# PART 1 -- Monte Carlo under the stated model
# ---------------------------------------------------------------------------
def part1(reps=4000, n=1000):
    se, f0 = 0.90, 0.0
    p1 = p1_of(0.60, P0)
    dt = delta_tip(0.60, P0, se, f0) / 100.0

    rows = []
    for mult in (0.5, 0.8, 1.0, 1.2, 1.5):
        d = dt * mult
        r0 = observed_risk_effective(P0, se, f0)
        r1 = observed_risk_effective(p1, se, f0 + d)
        empirical = []
        for _ in range(reps):
            x0 = sum(1 for _ in range(n) if random.random() < r0)
            x1 = sum(1 for _ in range(n) if random.random() < r1)
            rr0, rr1 = x0 / n, x1 / n
            if 0 < rr0 < 1 and 0 < rr1 < 1:
                empirical.append(or_from_risks(rr0, rr1))
        empirical.sort()
        rows.append({
            "delta_multiple_of_tip": mult,
            "delta_pp": 100 * d,
            "analytic_observed_or": or_from_risks(r0, r1),
            "mc_median_observed_or": empirical[len(empirical) // 2],
            "mc_ci2.5": empirical[int(0.025 * len(empirical))],
            "mc_ci97.5": empirical[int(0.975 * len(empirical))],
        })
    return {"model": "stated model (equal effective se, f vs f+Delta)",
            "se": se, "f0": f0, "p0": P0, "p1": p1, "true_or": 0.60,
            "delta_tip_pp": 100 * dt, "reps": reps, "n_per_arm": n, "rows": rows}


# ---------------------------------------------------------------------------
# PART 2 -- the boundary-error model, and the conversion to Delta_tip units
# ---------------------------------------------------------------------------
def part2():
    a0, b0 = 0.15, 0.10        # reference arm boundary error, illustrative
    se0, f0 = effective(GRADES_REF, a0, b0)
    se1_base, f1_base = effective(GRADES_IDX, a0, b0)

    # (i) what the closed form predicts using the reference arm's effective values
    dt_closed = delta_tip(TRUE_OR, P0, se0, f0)

    # (ii) sweep the index arm's upgrade rate a1, holding b1 = b0, and find the
    #      value at which the observed OR crosses 1
    def obs_or_at_a1(a1, b1=b0):
        r0 = observed_risk_grade(GRADES_REF, a0, b0)
        r1 = observed_risk_grade(GRADES_IDX, a1, b1)
        return or_from_risks(r0, r1)

    lo, hi = a0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if obs_or_at_a1(mid) < 1.0:
            lo = mid
        else:
            hi = mid
    a1_star = (lo + hi) / 2
    se1_star, f1_star = effective(GRADES_IDX, a1_star, b0)
    delta_needed = 100.0 * (f1_star - f0)

    # (iii) what a common (non-differential) process alone does
    r0_common = observed_risk_grade(GRADES_REF, a0, b0)
    r1_common = observed_risk_grade(GRADES_IDX, a0, b0)
    or_common = or_from_risks(r0_common, r1_common)
    true_or_from_dist = or_from_risks(P0, P1)

    # effective values under the common process, and the "Delta" they imply
    delta_apparent = 100.0 * (f1_base - f0)

    return {
        "reference_arm": {"a": a0, "b": b0, "se": se0, "f": f0,
                          "p_true": P0, "observed_risk": r0_common},
        "index_arm_common_process": {"a": a0, "b": b0, "se": se1_base, "f": f1_base,
                                     "p_true": P1, "observed_risk": r1_common},
        "true_or": true_or_from_dist,
        "or_under_common_process": or_common,
        "bias_from_common_process_pp_of_or": 100 * (or_common - true_or_from_dist),
        "delta_apparent_from_common_process_pp": delta_apparent,
        "delta_tip_closed_form_pp": dt_closed,
        "index_upgrade_rate_for_null": a1_star,
        "index_effective_se_at_null": se1_star,
        "index_effective_f_at_null": f1_star,
        "delta_needed_for_null_pp": delta_needed,
        "ratio_needed_to_closed_form": delta_needed / dt_closed if dt_closed else None,
        "ratio_apparent_to_closed_form": (delta_apparent / dt_closed) if dt_closed else None,
    }


# ---------------------------------------------------------------------------
# PART 3 -- two-way asymmetry: the tipping surface
# ---------------------------------------------------------------------------
def part3():
    a0, b0 = 0.15, 0.10
    r0 = observed_risk_grade(GRADES_REF, a0, b0)
    grid = []
    for a1_pp in range(5, 81, 5):
        for b1_pp in range(0, 41, 5):
            a1, b1 = a1_pp / 100.0, b1_pp / 100.0
            r1 = observed_risk_grade(GRADES_IDX, a1, b1)
            grid.append({"a1": a1, "b1": b1,
                         "observed_or": or_from_risks(r0, r1),
                         "delta_effective_pp": 100.0 * (effective(GRADES_IDX, a1, b1)[1] - effective(GRADES_REF, a0, b0)[1])})
    # the one-way slice the manuscript reports is b1 = b0
    slice_ = [g for g in grid if abs(g["b1"] - b0) < 1e-9]

    # The OR = 1 contour: for each b1, the a1 that yields a null. This shows how
    # much the required false-positive asymmetry moves once the OTHER direction
    # of error is allowed to differ too -- i.e. how wrong the one-way answer is
    # when the assumption it rests on fails.
    contour = []
    for b1_pp in range(0, 41, 5):
        b1 = b1_pp / 100.0

        def f(a1):
            return or_from_risks(r0, observed_risk_grade(GRADES_IDX, a1, b1))

        if f(0.0) >= 1.0 or f(1.0) < 1.0:
            contour.append({"b1": b1, "a1_for_null": None,
                            "delta_effective_pp": None})
            continue
        lo, hi = 0.0, 1.0
        for _ in range(200):
            mid = (lo + hi) / 2
            if f(mid) < 1.0:
                lo = mid
            else:
                hi = mid
        a1s = (lo + hi) / 2
        contour.append({
            "b1": b1,
            "a1_for_null": a1s,
            "delta_effective_pp": 100.0 * (effective(GRADES_IDX, a1s, b1)[1]
                                           - effective(GRADES_REF, a0, b0)[1]),
            "se_index": effective(GRADES_IDX, a1s, b1)[0],
        })
    return {"a0": a0, "b0": b0, "reference_observed_risk": r0,
            "grid": grid, "one_way_slice": slice_, "or1_contour": contour}


# ---------------------------------------------------------------------------
# PART 4 -- robustness: does the artefact depend on the grade-distribution shape?
# ---------------------------------------------------------------------------
def _dist(poor_risk, concentration, grades=(0, 1, 2, 3, 4, 5, 6)):
    """Build a 7-grade distribution with a given poor-outcome risk.

    `concentration` (0..1) controls how much of the favourable mass sits at
    mRS 2 rather than mRS 0-1. Varying it changes P(true=2 | true<=2), which is
    the conversion factor between boundary error and effective false positives.
    """
    fav = {0: (1 - concentration) * 0.4, 1: (1 - concentration) * 0.6,
           2: concentration}
    poor = {3: 0.28, 4: 0.28, 5: 0.20, 6: 0.24}
    d = {}
    for g in grades:
        if g <= 2:
            d[g] = (1 - poor_risk) * fav[g]
        else:
            d[g] = poor_risk * poor[g]
    return d


def part4():
    a0, b0 = 0.15, 0.10
    rows = []
    for conc_ref in (0.30, 0.40, 0.50, 0.60):
        for conc_idx in (0.30, 0.40, 0.50, 0.60):
            dref = _dist(0.55, conc_ref)
            didx = _dist(0.423, conc_idx)
            se0, f0 = effective(dref, a0, b0)
            se1, f1 = effective(didx, a0, b0)
            r0 = observed_risk_grade(dref, a0, b0)
            r1 = observed_risk_grade(didx, a0, b0)
            dt = delta_tip(TRUE_OR, P0, se0, f0)
            rows.append({
                "conc_ref": conc_ref, "conc_idx": conc_idx,
                "se_ref": se0, "f_ref": f0, "se_idx": se1, "f_idx": f1,
                "apparent_delta_pp": 100 * (f1 - f0),
                "apparent_delta_as_pct_of_tip": 100 * (f1 - f0) / (dt / 100.0) if dt else None,
                "observed_or_common_process": or_from_risks(r0, r1),
                "true_or": or_from_risks(
                    sum(v for g, v in dref.items() if g >= 3),
                    sum(v for g, v in didx.items() if g >= 3)),
                "delta_tip_closed_form_pp": dt,
            })
    ad = [abs(r["apparent_delta_pp"]) for r in rows]
    return {"rows": rows,
            "apparent_delta_abs_range_pp": [min(ad), max(ad)],
            "note": "apparent Delta is what an analyst would compute from the "
                    "arm-specific effective rates; it is non-zero although the "
                    "measurement process is identical in the two arms."}


def main():
    out = {"p0": P0, "p1": P1, "true_or": TRUE_OR,
           "grades_ref": GRADES_REF, "grades_idx": GRADES_IDX,
           "part1_stated_model": part1(),
           "part2_boundary_model": part2(),
           "part3_surface": part3(),
           "part4_shape_robustness": part4()}

    os.makedirs("reports", exist_ok=True)
    with io.open("reports/simulation_tipping_point.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    # ---- console report ----
    print("=" * 74)
    print("SETUP")
    print("=" * 74)
    print(f"  reference arm true poor-outcome risk p0 = {P0:.3f}")
    print(f"  index arm     true poor-outcome risk p1 = {P1:.3f}")
    print(f"  true odds ratio = {TRUE_OR:.4f}")

    p1_ = out["part1_stated_model"]
    print()
    print("=" * 74)
    print("PART 1 -- Monte Carlo under the stated model (numerical check of algebra)")
    print("=" * 74)
    print(f"  se = {p1_['se']}, f0 = {p1_['f0']}, Delta_tip = {p1_['delta_tip_pp']:.2f} pp")
    print(f"  {'mult':>5} {'Delta(pp)':>10} {'analytic OR':>12} {'MC median':>11} {'95% MC CI':>22}")
    for r in p1_["rows"]:
        print(f"  {r['delta_multiple_of_tip']:>5.1f} {r['delta_pp']:>10.2f} "
              f"{r['analytic_observed_or']:>12.4f} {r['mc_median_observed_or']:>11.4f} "
              f"{'[' + format(r['mc_ci2.5'], '.3f') + ', ' + format(r['mc_ci97.5'], '.3f') + ']':>22}")

    p2 = out["part2_boundary_model"]
    print()
    print("=" * 74)
    print("PART 2 -- boundary-error model (error where it actually occurs)")
    print("=" * 74)
    print(f"  reference arm  a={p2['reference_arm']['a']}, b={p2['reference_arm']['b']}"
          f"  ->  se={p2['reference_arm']['se']:.4f}, f={p2['reference_arm']['f']:.4f}")
    ic = p2["index_arm_common_process"]
    print(f"  index arm, SAME process  a={ic['a']}, b={ic['b']}"
          f"  ->  se={ic['se']:.4f}, f={ic['f']:.4f}   <-- differs anyway")
    print()
    print(f"  true OR                                 = {p2['true_or']:.4f}")
    print(f"  OR under the common (non-differential) process = {p2['or_under_common_process']:.4f}")
    print(f"  apparent 'Delta' from that common process = {p2['delta_apparent_from_common_process_pp']:.2f} pp")
    print()
    print(f"  closed form Delta_tip (reference effective se, f) = {p2['delta_tip_closed_form_pp']:.2f} pp")
    print(f"  upgrade rate in index arm needed for a null      = {p2['index_upgrade_rate_for_null']:.4f}")
    print(f"  effective Delta at that point                    = {p2['delta_needed_for_null_pp']:.2f} pp")
    print(f"  ratio needed / closed form                       = {p2['ratio_needed_to_closed_form']:.3f}")

    print()
    print("=" * 74)
    print("PART 3 -- two-way asymmetry (one-way slice is the manuscript's)")
    print("=" * 74)
    print(f"  {'a1':>6} {'b1':>6} {'observed OR':>12} {'Delta_eff(pp)':>14}")
    for g in out["part3_surface"]["one_way_slice"]:
        print(f"  {g['a1']:>6.2f} {g['b1']:>6.2f} {g['observed_or']:>12.4f} {g['delta_effective_pp']:>14.2f}")

    print()
    print("=" * 74)
    print("PART 3b -- OR = 1 contour: what happens when the OTHER direction of")
    print("           error is also allowed to differ (one-way assumption fails)")
    print("=" * 74)
    print(f"  {'b1':>6} {'a1 for null':>13} {'Delta_eff(pp)':>15} {'se_index':>10}")
    for c in out["part3_surface"]["or1_contour"]:
        if c["a1_for_null"] is None:
            print(f"  {c['b1']:>6.2f} {'no null reachable':>13}")
        else:
            print(f"  {c['b1']:>6.2f} {c['a1_for_null']:>13.4f} "
                  f"{c['delta_effective_pp']:>15.2f} {c['se_index']:>10.4f}")

    p4 = out["part4_shape_robustness"]
    print()
    print("=" * 74)
    print("PART 4 -- robustness of the artefact to grade-distribution shape")
    print("=" * 74)
    print(f"  {'conc_ref':>9} {'conc_idx':>9} {'apparent Delta(pp)':>19} "
          f"{'as % of tip':>12} {'OR (common proc.)':>18}")
    for r in p4["rows"]:
        print(f"  {r['conc_ref']:>9.2f} {r['conc_idx']:>9.2f} "
              f"{r['apparent_delta_pp']:>19.2f} "
              f"{r['apparent_delta_as_pct_of_tip']:>11.1f}% "
              f"{r['observed_or_common_process']:>18.4f}")
    print(f"\n  |apparent Delta| across shapes: "
          f"{p4['apparent_delta_abs_range_pp'][0]:.2f} to "
          f"{p4['apparent_delta_abs_range_pp'][1]:.2f} pp")

    print()
    print("saved -> reports/simulation_tipping_point.json")


if __name__ == "__main__":
    main()
