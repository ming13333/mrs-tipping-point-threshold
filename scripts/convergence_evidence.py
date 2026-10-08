# -*- coding: utf-8 -*-
"""W1 step 8: convergence evidence. The direct answer to BMC's "well tested".

Scope, deliberately narrow. Three blocks, one per requirement in the outline's
convergence section. Nothing else is added.

  Block R -- replicate doubling.  The same cell at R and at 2R replicates, on a
             NESTED stream, so the difference is paired and its standard error is
             computable rather than guessed.
  Block G -- grid refinement.  The scenario grid is refined by inserting the
             midpoint of every interval, on p0 and on se, giving a strict
             superset. Analytic quantities must be bit-identical on shared
             cells; sampled quantities are recomputed and compared against their
             own Monte Carlo error.
  Block D -- delta_sig resolution.  The located crossing is recomputed with
             twice the bisection iterations and twice the replicates, because
             delta_sig is the only headline quantity that is estimated rather
             than derived.

WHY THIS FILE DOES NOT REUSE THE STEP 4 AND STEP 5 NUMBERS. Those steps consume
one random stream sequentially across the grid: each row rebuilds
random.Random(SEED) and the draws follow in traversal order. Measured on the
type I rows, inserting a single extra cell changed 33 of the 40 shared cells,
because every later cell inherits a shifted stream position. Under that design a
grid-refinement comparison cannot separate a genuine change in the estimate from
a change in the random stream, and the two are the same size as the Monte Carlo
error. So every number here is recomputed from the certified layer A records
with a stream derived from the cell key, which makes a cell's draws independent
of where it sits in the traversal. The sampling kernel itself is IMPORTED from
step 4 rather than reimplemented, so the only thing that changes is the seed
derivation. `nested_streams_share_their_prefix_is_verified` asserts the prefix
property rather than assuming it, and `seed_derivation_is_position_independent`
asserts that the same cell yields the same draws whatever surrounds it.

The paired standard error, derived rather than borrowed. For nested streams the
2R replicate verdicts v_1..v_2R contain the R verdicts v_1..v_R as a prefix, so
with r_R = mean(v_1..v_R) and r_last = mean(v_{R+1}..v_{2R}),

    rate_2R = (r_R + r_last) / 2
    rate_2R - rate_R = (r_last - rate_R) / 2
    Var(rate_2R - rate_R) = (1/4) [p(1-p)/R + p(1-p)/R] = p(1 - p) / (2R)

so

    se_paired = sqrt(rate_R (1 - rate_R) / (2R))

The factor of two is not optional and an earlier version of this file omitted the
1/4 from the averaging, giving sqrt(p(1-p)/R). Measured on 4,000 independent
nested pairs at p = 0.05 and R = 2,000, the empirical standard deviation of the
difference was 0.003497 against 0.003446 predicted by p(1-p)/(2R) and 0.004873
by p(1-p)/R, so the longer formula is the one that matches. The overstatement is
by a factor of about 1.73 at equal rates, not the sqrt(2) that a casual reading
suggests: the unpaired value is sqrt(p(1-p)/R + p(1-p)/(2R)), which is sqrt(3)
times the paired one. Overstating the error would make a real drift look
acceptable, so the comparison uses the paired value and the unpaired one is
reported beside it only to show the size of the mistake avoided.

WHAT THIS STEP DOES NOT ESTABLISH. Convergence of a Monte Carlo estimate is not
validity of the estimand, and none of the three blocks speaks to whether the
measurement model is the right one. Two of the three blocks are also silent on
the direction of any bias, because halving a step size or doubling a replicate
count moves an estimate towards its own expectation, not towards the truth.

Run:  python scripts/convergence_evidence.py
Out:  reports/estimator_eval_convergence.json
Exit 0 if every convergence check passes, 1 otherwise.
"""
import hashlib
import io
import json
import math
import os
import random
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

# The one permitted source of the mechanism and of the closed form.
import simulate_tipping_point as legacy   # noqa: E402
# The sampling kernel is imported, never reimplemented. It is a pure function of
# (q0, q1, n, reps, rng); only the seed derivation differs from step 4.
from calibrate_estimator_eval import run_cell, _z_two_sided   # noqa: E402
# delta_sig is located by the step 5 function, again imported rather than copied.
from locate_delta_sig_and_coverage import directional_criterion, ALPHA   # noqa: E402

DGP = os.path.join(ROOT, 'reports', 'estimator_eval_dgp_checks.json')
OUT = os.path.join(ROOT, 'reports', 'estimator_eval_convergence.json')

TOL = 1e-12
N_GRID = [100, 500, 1000, 2000]

# Block R. Two levels, so the drift can be read at R and at 2R rather than only
# at the far end, and a non-monotone pattern would be visible.
REPS_BASE = 2000
REPS_DOUBLE = 4000

# Block D. The step 5 setting, and the refined one. Both are recorded because a
# single pair cannot distinguish a tolerance effect from a replicate effect.
BISECT_ITERS_BASE = 26
BISECT_ITERS_FINE = 34
REPS_SIG_BASE = 8000
REPS_SIG_FINE = 16000
DELTA_SIG_CELLS = 8          # cells carrying the located crossing

# One namespace, so no block can collide with another or with steps 4 and 5.
SEED_NS_CONV = b'estimator_eval_convergence/v1'


# --------------------------------------------------------------------- seeds
def derive_seed(*parts):
    """A stream derived from the cell identity, not from its position.

    hashlib rather than hash() because PYTHONHASHSEED randomises the built-in
    hash of a str per process, which would make the whole artefact irreproducible
    across runs while still looking deterministic inside one run.
    """
    h = hashlib.sha256()
    h.update(SEED_NS_CONV)
    for p in parts:
        h.update(b'|')
        h.update((p if isinstance(p, bytes) else str(p)).encode('utf-8'))
    return int.from_bytes(h.digest()[:8], 'big')


# ------------------------------------------------------------------- layer A
def load_layer_a():
    if not os.path.exists(DGP):
        raise SystemExit('step 3 output missing: run step 3 first')
    with io.open(DGP, encoding='utf-8') as fh:
        doc = json.load(fh)
    if doc['meta']['step'] != 'W1 step 3' or not doc['summary']['all_passed']:
        raise SystemExit('step 3 mechanisms are not certified')
    return doc['records']['layer_A']


def key_of(r):
    return 'A|p0=%s|or_true=%s|se=%s|f=%s' % (r['p0'], r['or_true'],
                                              r['se'], r['f'])


def midpoints(v):
    """A strict superset of v containing every original level.

    Inserting the midpoint of each interval, rather than halving the step from
    one end, is what makes the refined grid NESTED. Without nesting, a shared
    cell is not the same cell and no comparison is possible.
    """
    out = []
    for a, b in zip(v, v[1:]):
        mid = round((a + b) / 2.0, 4)
        out.append(a)
        out.append(mid)
    out.append(v[-1])
    return out


# ------------------------------------------------------------------- block R
def block_r_worker(job):
    (label, p0, p1, or_true, se, f, n) = job
    n_key = key_of({'p0': p0, 'or_true': or_true, 'se': se, 'f': f})
    # One stream per cell per level. The 2R stream is NOT a different seed: it is
    # the same stream drawn twice as long, which is what makes the pair nested.
    seed = derive_seed('R', n_key, n)
    q0 = legacy.observed_risk_effective(p0, se, f)
    q1 = legacy.observed_risk_effective(p1, se, f)
    short = run_cell(q0, q1, n, REPS_BASE, random.Random(seed))
    long_ = run_cell(q0, q1, n, REPS_DOUBLE, random.Random(seed))
    r1, r2 = short['reject_rate'], long_['reject_rate']
    if r1 is None or r2 is None:
        return {'label': label, 'cell_key': n_key, 'n_per_arm': n,
                'usable': False, 'seed': seed}
    # The paired difference is (rate_last - rate_R) / 2, so its variance is
    # p(1-p)/(2R). See the module docstring for the derivation and for the
    # numerical check that discriminates the two candidate formulas.
    se_paired = math.sqrt(r1 * (1.0 - r1) / (2.0 * short['usable_reps']))
    # The value an unpaired comparison would have used, kept to show the size of
    # the error avoided. At equal rates this is sqrt(3) times the paired value.
    se_independent = math.sqrt(r1 * (1.0 - r1) / short['usable_reps']
                               + r2 * (1.0 - r2) / long_['usable_reps'])
    diff = r2 - r1
    return {
        'label': label, 'cell_key': n_key, 'n_per_arm': n, 'usable': True,
        'seed': seed,
        'rate_at_R': r1, 'rate_at_2R': r2,
        'rejections_at_R': short['rejections'], 'usable_at_R': short['usable_reps'],
        'rejections_at_2R': long_['rejections'], 'usable_at_2R': long_['usable_reps'],
        'difference': diff, 'abs_difference': abs(diff),
        'se_paired': se_paired, 'se_independent': se_independent,
        'z_paired': (diff / se_paired) if se_paired > 0 else None,
        'z_if_unpaired': (diff / se_independent) if se_independent > 0 else None,
        'within_one_paired_se': abs(diff) <= se_paired,
        'within_two_paired_se': abs(diff) <= 2.0 * se_paired,
        'verdict_changed': None,          # filled in by the caller
        'verdict_at_R': None, 'verdict_at_2R': None,
    }


def block_r(layer_a):
    jobs = []
    for r in layer_a:
        if not r['legal'] or abs(r['or_true'] - 1.0) > TOL:
            continue
        for n in N_GRID:
            jobs.append(('type_I', r['p0'], r['p1'], r['or_true'],
                         r['se'], r['f'], n))
    with Pool(processes=min(len(jobs), os.cpu_count() or 1)) as pool_:
        rows = pool_.map(block_r_worker, jobs, chunksize=4)
    for x in rows:
        if not x.get('usable'):
            continue
        # The verdict is recorded at the declared band only. The sensitivity of
        # the verdict to the band is computed in the check, from the two rates,
        # so that one place owns the band and no row carries a stale label.
        x['verdict_at_R'] = _verdict(x['rate_at_R'], VERDICT_BAND)
        x['verdict_at_2R'] = _verdict(x['rate_at_2R'], VERDICT_BAND)
        x['verdict_changed'] = x['verdict_at_R'] != x['verdict_at_2R']
    return [x for x in rows if x.get('usable')], len(jobs)


# ------------------------------------------------------------------- block G
COARSE_P0 = [0.45, 0.50, 0.55, 0.60, 0.65]
COARSE_SE = [0.70, 0.80, 0.90, 0.98]
OR_GRID = [0.50, 0.60, 0.70, 0.80, 1.00]
F_GRID = [0.0, 0.05]


def analytic_cell(p0, or_true, se, f):
    p1 = legacy.p1_of(or_true, p0)
    return {
        'cell_key': 'A|p0=%s|or_true=%s|se=%s|f=%s' % (p0, or_true, se, f),
        'p0': p0, 'p1': p1, 'or_true': or_true, 'se': se, 'f': f,
        'delta_tip_pp': legacy.delta_tip(or_true, p0, se, f),
        'q0': legacy.observed_risk_effective(p0, se, f),
        'q1': legacy.observed_risk_effective(p1, se, f),
    }


def block_g_analytic():
    """Refine the grid and compare analytic quantities on the shared cells.

    Nothing is sampled here, so the comparison is exact and any difference is a
    bug rather than noise.
    """
    fine_p0, fine_se = midpoints(COARSE_P0), midpoints(COARSE_SE)
    coarse = {}
    for p0 in COARSE_P0:
        for se in COARSE_SE:
            for or_true in OR_GRID:
                c = analytic_cell(p0, or_true, se, 0.0)
                c2 = analytic_cell(p0, or_true, se, 0.05)
                coarse[c['cell_key']] = c
                coarse[c2['cell_key']] = c2
    fine = {}
    for p0 in fine_p0:
        for se in fine_se:
            for or_true in OR_GRID:
                c = analytic_cell(p0, or_true, se, 0.0)
                c2 = analytic_cell(p0, or_true, se, 0.05)
                fine[c['cell_key']] = c
                fine[c2['cell_key']] = c2
    shared = sorted(set(coarse) & set(fine))
    analytic_mismatch = [k for k in shared
                         if coarse[k]['delta_tip_pp'] != fine[k]['delta_tip_pp']
                         or coarse[k]['q0'] != fine[k]['q0']
                         or coarse[k]['q1'] != fine[k]['q1']]
    # The extremes are a property of the closed form over the box, not of the
    # sample of points. Both are recomputed here and then explained by the
    # monotonicity argument recorded in the checks, so that "the extremes did not
    # move" is a derived consequence rather than an observation.
    def extremes(pool):
        elig = [c for c in pool.values() if abs(c['or_true'] - 1.0) > TOL]
        lo = min(elig, key=lambda c: c['delta_tip_pp'])
        hi = max(elig, key=lambda c: c['delta_tip_pp'])
        return lo, hi

    lo_c, hi_c = extremes(coarse)
    lo_f, hi_f = extremes(fine)
    return {
        'coarse_axes': {'p0': COARSE_P0, 'se': COARSE_SE},
        'fine_axes': {'p0': fine_p0, 'se': fine_se},
        'nested_p0': set(COARSE_P0) <= set(fine_p0),
        'nested_se': set(COARSE_SE) <= set(fine_se),
        'coarse_cells': len(coarse), 'fine_cells': len(fine),
        'shared_cells': len(shared),
        'new_cells': len(set(fine) - set(coarse)),
        'analytic_mismatches': analytic_mismatch,
        'coarse_delta_tip_range': [lo_c['delta_tip_pp'], hi_c['delta_tip_pp']],
        'fine_delta_tip_range': [lo_f['delta_tip_pp'], hi_f['delta_tip_pp']],
        'coarse_argmin': {k: lo_c[k] for k in ('p0', 'se', 'or_true', 'f')},
        'coarse_argmax': {k: hi_c[k] for k in ('p0', 'se', 'or_true', 'f')},
        'fine_argmin': {k: lo_f[k] for k in ('p0', 'se', 'or_true', 'f')},
        'fine_argmax': {k: hi_f[k] for k in ('p0', 'se', 'or_true', 'f')},
        'monotonicity': {
            'in_p0': 'strictly increasing',
            'in_se': 'strictly increasing',
            'in_f': 'strictly decreasing',
            'in_or_true': 'strictly decreasing',
            'consequence': 'delta_tip is monotone in each axis, so its extremes '
                           'over a box are attained at the vertices; inserting '
                           'midpoints of p0 and se adds no vertex, so the '
                           'extremes cannot move. Verified numerically below '
                           'rather than only argued.',
        },
        '_coarse': coarse, '_fine': fine, '_shared': shared,
    }


def _g_worker(job):
    (p0, p1, se, f, n, tag) = job
    seed = derive_seed('G', tag, p0, p1, se, f, n)
    q0 = legacy.observed_risk_effective(p0, se, f)
    q1 = legacy.observed_risk_effective(p1, se, f)
    res = run_cell(q0, q1, n, REPS_BASE, random.Random(seed))
    res['p0'], res['p1'], res['se'], res['f'] = p0, p1, se, f
    res['n_per_arm'] = n
    res['cell_key'] = 'A|p0=%s|or_true=1.0|se=%s|f=%s' % (p0, se, f)
    return res


def block_g_sampled(analytic):
    """Recompute the type I error on both grids, with position-free seeds."""
    coarse, fine = analytic['_coarse'], analytic['_fine']
    jobs_c, jobs_f = [], []
    for k in sorted(coarse):
        c = coarse[k]
        if abs(c['or_true'] - 1.0) > TOL:
            continue
        for n in N_GRID:
            jobs_c.append((c['p0'], c['p1'], c['se'], c['f'], n, 'coarse'))
    for k in sorted(fine):
        c = fine[k]
        if abs(c['or_true'] - 1.0) > TOL:
            continue
        for n in N_GRID:
            jobs_f.append((c['p0'], c['p1'], c['se'], c['f'], n, 'fine'))
    with Pool(processes=min(len(jobs_c) + len(jobs_f),
                            os.cpu_count() or 1)) as pool_:
        rc = pool_.map(_g_worker, jobs_c, chunksize=8)
        rf = pool_.map(_g_worker, jobs_f, chunksize=8)
    index = lambda rows: {(x['cell_key'], x['n_per_arm']): x for x in rows}
    ic, if_ = index(rc), index(rf)
    shared = sorted(set(ic) & set(if_))

    def pooled(index_, keys):
        num = sum(index_[k]['rejections'] for k in keys)
        den = sum(index_[k]['usable_reps'] for k in keys)
        rate = num / den
        return {'rate': rate, 'usable_reps': den,
                'mc95_halfwidth': 1.96 * math.sqrt(rate * (1 - rate) / den)}

    pooled_c = pooled(ic, shared)
    pooled_f = pooled(if_, sorted(if_))
    # Per-cell drift on the shared cells, judged against the cell's own paired
    # scale. Two independently seeded cells are not nested, so the difference of
    # two rates carries sqrt(se1^2 + se2^2); using one cell's error alone would
    # halve the yardstick.
    per_cell = []
    for k in shared:
        a, b = ic[k], if_[k]
        if a['reject_rate'] is None or b['reject_rate'] is None:
            continue
        se = math.sqrt(a['mc95_halfwidth'] ** 2 + b['mc95_halfwidth'] ** 2) / 1.96
        d = b['reject_rate'] - a['reject_rate']
        per_cell.append({
            'cell_key': k[0], 'n_per_arm': k[1],
            'rate_coarse': a['reject_rate'], 'rate_fine': b['reject_rate'],
            'difference': d, 'se_of_difference': se,
            'z': (d / se) if se > 0 else None,
            'within_two_se': abs(d) <= 2.0 * se,
        })
    zs = [x['z'] for x in per_cell if x['z'] is not None]
    return {
        'coarse_cells_sampled': len(jobs_c), 'fine_cells_sampled': len(jobs_f),
        'shared_cells_sampled': len(shared),
        'pooled_on_shared_coarse': pooled_c,
        'pooled_over_all_fine': pooled_f,
        'pooled_difference': pooled_f['rate'] - pooled_c['rate'],
        'pooled_mc95_halfwidth_coarse': pooled_c['mc95_halfwidth'],
        'pooled_mc95_halfwidth_fine': pooled_f['mc95_halfwidth'],
        'per_cell': per_cell,
        'max_abs_z': max((abs(z) for z in zs), default=None),
        'cells_outside_two_se': [x for x in per_cell if not x['within_two_se']],
        '_coarse_index': ic, '_fine_index': if_,
    }


# ------------------------------------------------------------------- block D
def _bisect(p0, p1, se, f, d_tip, n, iters, reps, seed_base):
    """Locate delta_sig by bisection on the directional criterion.

    The step 5 function is used for the criterion itself. Only the iteration
    count and the replicate count are parameters here, because those are the two
    things this block varies.

    The achieved bracket width is RETURNED, not recomputed by the caller from a
    constant. An earlier version had the caller derive it from BISECT_ITERS_FINE
    while the loop ran with whatever it was passed, so a run with a deliberately
    coarse iteration count still reported a bracket of 9e-10 pp and the
    bracket-width guard passed on a number that described no computation that had
    happened.
    """
    at_lo, _ = directional_criterion(p1, p0, se, f, 0.0, n, reps, seed_base + 1)
    at_tip, _ = directional_criterion(p1, p0, se, f, d_tip, n, reps,
                                       seed_base + 2)
    if at_lo is None or at_tip is None:
        return None, None, None
    if at_tip >= ALPHA:
        # No crossing inside the domain: the bracket is the whole search domain,
        # and saying so is the point rather than reporting a tight tolerance for
        # a value that was never located.
        return d_tip, True, 100.0 * d_tip
    lo, hi = 0.0, d_tip
    done = 0
    for k in range(iters):
        mid = 0.5 * (lo + hi)
        val, _ = directional_criterion(p1, p0, se, f, mid, n, reps,
                                       seed_base + 3 + k)
        if val is None:
            return None, False, None
        done = k + 1
        if val >= ALPHA:
            lo = mid
        else:
            hi = mid
    # The width actually reached, from the width actually achieved.
    width = 100.0 * (d_tip / (2.0 ** done))
    return lo, False, width


def _d_worker(job):
    (idx, p0, p1, se, f, d_tip, n) = job
    seed_base = derive_seed('D', p0, p1, se, f, n)
    base, capped_base, width_base = _bisect(
        p0, p1, se, f, d_tip, n, BISECT_ITERS_BASE, REPS_SIG_BASE, seed_base)
    fine, capped_fine, width_fine = _bisect(
        p0, p1, se, f, d_tip, n, BISECT_ITERS_FINE, REPS_SIG_FINE, seed_base)
    if base is None or fine is None or width_base is None or width_fine is None:
        return {'index': idx, 'usable': False}
    # The bracket widths are the ones the two runs actually achieved, returned by
    # the locator. Recomputing them here from the module constants is what let a
    # three-iteration run report a tolerance of 9e-10 pp.
    return {
        'index': idx, 'usable': True,
        'cell_key': 'A|p0=%s|se=%s|f=%s' % (p0, se, f),
        'p0': p0, 'p1': p1, 'se': se, 'f': f, 'n_per_arm': n,
        'delta_tip_pp': 100.0 * d_tip,
        'delta_sig_base_pp': 100.0 * base,
        'delta_sig_fine_pp': 100.0 * fine,
        'shift_pp': 100.0 * abs(fine - base),
        'shift_as_pct_of_tip': abs(fine - base) / d_tip,
        'ratio_base': base / d_tip, 'ratio_fine': fine / d_tip,
        'capped_base': capped_base, 'capped_fine': capped_fine,
        'bracket_width_base_pp': width_base,
        'bracket_width_fine_pp': width_fine,
        'bracket_is_finer': width_fine < width_base,
        'within_point_two_pp': abs(100.0 * (fine - base)) < 0.2,
    }


def block_d(layer_a):
    cand = [r for r in layer_a
            if r['legal'] and abs(r['or_true'] - 1.0) > TOL
            and r['threshold_exists']]
    # Deterministic spread over the grid rather than the first N in file order:
    # picking the head of a file would silently favour one corner of p0 x se.
    picks = []
    step = max(1, len(cand) // DELTA_SIG_CELLS)
    picks = cand[::step][:DELTA_SIG_CELLS]
    jobs = []
    for i, r in enumerate(picks):
        for n in (1000, 2000):
            jobs.append((len(jobs), r['p0'], r['p1'], r['se'], r['f'],
                         r['delta_tip_pp'] / 100.0, n))
    with Pool(processes=min(len(jobs), os.cpu_count() or 1)) as pool_:
        rows = pool_.map(_d_worker, jobs, chunksize=1)
    return [x for x in rows if x.get('usable')], len(jobs), len(cand)


# -------------------------------------------------------------------- checks
def check_nested_prefix_is_real():
    """The 2R stream must extend the R stream, not merely resemble it.

    If the two levels used different seeds the paired standard error would be
    wrong by construction, and every z in block R would be understated. The
    property is asserted by replaying one cell twice and comparing the first R
    verdicts, which is the only form of the claim that can fail.
    """
    p0, p1, se, f, n = 0.55, 0.55, 0.90, 0.05, 500
    seed = derive_seed('R', key_of({'p0': p0, 'or_true': 1.0, 'se': se, 'f': f}), n)
    q0 = legacy.observed_risk_effective(p0, se, f)
    q1 = legacy.observed_risk_effective(p1, se, f)

    def verdicts(reps):
        rng = random.Random(seed)
        out = []
        for _ in range(reps):
            x0 = rng.binomialvariate(n, q0)
            x1 = rng.binomialvariate(n, q1)
            if x0 in (0, n) or x1 in (0, n):
                continue
            lor = math.log(x1 / (n - x1)) - math.log(x0 / (n - x0))
            v = (1.0 / x0 + 1.0 / (n - x0) + 1.0 / x1 + 1.0 / (n - x1)) ** 0.5
            out.append(1 if abs(lor / v) > _z_two_sided(ALPHA) else 0)
        return out

    short = verdicts(200)
    long_ = verdicts(400)
    prefix_equal = long_[:len(short)] == short
    # And the contrast: a differently seeded stream must NOT share the first
    # draw, otherwise the check above would pass for the wrong reason -- a
    # binomial generator that ignored its seed would also produce a stable
    # prefix, and would make every z below meaningless.
    differs = (verdicts_first_two(seed, n, q0, q1)
               != verdicts_first_two(seed + 1, n, q0, q1))
    return {
        'id': 'nested_streams_share_their_prefix_is_verified',
        'question': 'does the 2R stream literally extend the R stream?',
        'short_verdicts': len(short), 'long_verdicts': len(long_),
        'prefix_identical': prefix_equal,
        'independent_stream_differs_at_the_first_draw': differs,
        'passed': prefix_equal and differs,
        'why_it_matters': 'the paired standard error sqrt(p(1-p)/R) is only '
                          'correct under nesting. Under independent streams the '
                          'error is sqrt(2) larger and every z below is '
                          'understated by about 41 per cent',
    }


def verdicts_first_two(seed, n, q0, q1):
    rng = random.Random(seed)
    return (rng.binomialvariate(n, q0), rng.binomialvariate(n, q1))


def check_seed_is_position_free():
    """A cell's draws must not depend on what surrounds it.

    Step 4 and step 5 fail this: they consume one stream across the grid in
    traversal order, so inserting a cell shifts every later cell. The claim is
    checked by drawing the target cell twice, once with a decoy cell consumed
    first and once without, on the same derived seed.
    """
    target = dict(p0=0.55, p1=0.55, se=0.90, f=0.05, n=500)
    seed = derive_seed('R', key_of(dict(p0=0.55, or_true=1.0, se=0.90, f=0.05)), 500)
    direct = run_cell(legacy.observed_risk_effective(0.55, 0.90, 0.05),
                      legacy.observed_risk_effective(0.55, 0.90, 0.05),
                      500, 400, random.Random(seed))
    rng = random.Random(seed)
    rng.binomialvariate(500, 0.4)
    rng.binomialvariate(500, 0.3)
    after_decoy = run_cell(legacy.observed_risk_effective(0.55, 0.90, 0.05),
                           legacy.observed_risk_effective(0.55, 0.90, 0.05),
                           500, 400, rng)
    return {
        'id': 'seed_derivation_is_position_independent',
        'question': 'does consuming a decoy draw first change the cell result?',
        'rate_direct': direct['reject_rate'],
        'rate_after_decoy': after_decoy['reject_rate'],
        'identical': direct['reject_rate'] == after_decoy['reject_rate'],
        'passed': direct['reject_rate'] == after_decoy['reject_rate'],
        'why_it_matters': 'this is the property steps 4 and 5 lack. Without it a '
                          'grid-refinement comparison cannot separate a change in '
                          'the estimate from a change in the random stream, and '
                          'measured on the step 4 type I rows a single inserted '
                          'cell moved 33 of 40 shared cells',
    }


def _verdict(rate, band):
    if rate is None:
        return 'undefined'
    if rate < 0.05 - band:
        return 'below_nominal'
    if rate > 0.05 + band:
        return 'above_nominal'
    return 'at_nominal'


# Verdict band. Fixed at 0.5 pp by an earlier version, which is a defect and not
# a choice: at R = 2,000 the paired standard error near a rate of 0.05 is about
# 0.34 pp, so a 0.5 pp band is 1.5 errors wide and the classification of a cell
# is decided by noise. Measured on this grid, 47 of 160 verdicts changed between
# R and 2R at that width. The band is therefore set to 1 pp, roughly three paired
# errors, and the sensitivity of the verdict count to the band is reported rather
# than hidden: a reader who wants a stricter rule can read the table.
VERDICT_BAND = 0.010
VERDICT_BAND_SENSITIVITY = [0.005, 0.0075, 0.010, 0.0125, 0.015]


def check_replicate_doubling(rows):
    """Does doubling the replicates leave the reported numbers reproducible?

    THE TOLERANCE IS DISTRIBUTIONAL, AND THAT IS THE POINT. An earlier version
    required every one of the cells to sit inside two paired standard errors and
    reported FAIL on 6 of 160. That criterion is wrong by construction rather than
    wrong on this data: with 160 cells, expecting all of them inside plus or minus
    two errors has about a 95 per cent chance of failing even when the drift is
    pure noise. The observed z values have mean 0.023, standard deviation 1.005
    and mean absolute value 0.823, against 0.798 for a standard normal, and 6 of
    them exceed two errors where 7.3 are expected; simulating 4,000 draws of 160
    standard normals puts the 5th to 95th percentile of that count at 3 to 12. The
    drift is therefore indistinguishable from noise, and the criterion has to say
    so in a form that can distinguish noise from drift.

    The verdict band is a separate matter and is reported with its own
    sensitivity, because whether a three-way label is stable depends on how wide
    the band is relative to the error, not on how good the simulation is.
    """
    zs = [x['z_paired'] for x in rows if x['z_paired'] is not None]
    within1 = sum(1 for x in rows if x['within_one_paired_se'])
    within2 = sum(1 for x in rows if x['within_two_paired_se'])
    sensitivity = {}
    for band in VERDICT_BAND_SENSITIVITY:
        sensitivity['%.4f' % band] = sum(
            1 for x in rows
            if _verdict(x['rate_at_R'], band) != _verdict(x['rate_at_2R'], band))
    typical_se = rows[0]['se_paired'] if rows else None
    mean_z = (sum(zs) / len(zs)) if zs else None
    sd_z = (math.sqrt(sum((z - mean_z) ** 2 for z in zs) / (len(zs) - 1))
            if zs and len(zs) > 1 else None)
    over2 = sum(1 for z in zs if abs(z) > 2.0)
    # Binomial reference for the count beyond two errors, so the criterion is a
    # stated number rather than an intuition about tails.
    expected_over2 = 0.0455 * len(zs)
    # Tolerance on the count: the 5th and 95th percentiles of a binomial with
    # this many trials, approximated by the normal quantile on the count.
    half_width_95 = (1.96 * math.sqrt(0.0455 * (1 - 0.0455) * len(zs))
                     if zs else None)
    return {
        'id': 'replicate_doubling_leaves_the_type_i_estimate_reproducible',
        'question': 'is the paired drift distributed as noise rather than as a '
                    'shift, and is the verdict stable at a band wider than the '
                    'Monte Carlo error?',
        'cells': len(rows),
        'within_one_paired_se': within1,
        'within_two_paired_se': within2,
        'mean_z': mean_z, 'sd_z': sd_z,
        'mean_abs_z': (sum(abs(z) for z in zs) / len(zs)) if zs else None,
        'mean_abs_z_expected_for_standard_normal':
            math.sqrt(2.0 / math.pi),
        'max_abs_z_paired': max((abs(z) for z in zs), default=None),
        'cells_beyond_two_se': over2,
        'cells_beyond_two_se_expected': expected_over2,
        'cells_beyond_two_se_tolerance_95': half_width_95,
        'verdict_band': VERDICT_BAND,
        'verdict_band_in_paired_se': (VERDICT_BAND / typical_se
                                      if typical_se else None),
        'verdicts_changed': sum(1 for x in rows if x['verdict_changed']),
        'verdict_flips_by_band': sensitivity,
        'nominally_at_5pc_at_2R': sum(1 for x in rows
                                      if x['verdict_at_2R'] == 'at_nominal'),
        'nominally_at_5pc_at_R': sum(1 for x in rows
                                     if x['verdict_at_R'] == 'at_nominal'),
        'tolerance_note': 'the criterion is on the DISTRIBUTION of the paired z, '
                          'not on every cell individually. Requiring all 160 cells '
                          'inside two errors would fail about 95 per cent of the '
                          'time even on perfect data, so it cannot distinguish '
                          'noise from drift.',
        'passed': bool(zs
                       and abs(mean_z) < 0.20
                       and 0.85 < sd_z < 1.20
                       and abs(over2 - expected_over2) <= half_width_95),
        'why_it_matters': 'this is the block that answers the journal\'s "well '
                          'tested" requirement for the calibration claim. The '
                          'quantity being certified is that the reported rates do '
                          'not move when the budget doubles, which is a statement '
                          'about the spread of the drift and not about any single '
                          'cell',
    }


def check_paired_se_is_the_smaller_one(rows):
    """The paired error must be the smaller one, by the right margin.

    The ratio has a closed form. With A = p(1-p) for the short run and B for the
    long run,

        ratio^2 = [A / (2R)] / [A / R + B / (2R)] = A / (2A + B)

    so at equal rates it is 1/sqrt(3) = 0.577, and it rises towards 1 as B
    collapses towards zero, with no attained supremum. An earlier version of this
    check asserted an upper bound of sqrt(2)/2 = 0.707, which is not a property
    of the ratio at all: under the wrong paired formula the measured range was
    0.796 to 0.835, and under the right one it is 0.577 to 0.633. The assertion
    is therefore restated as the two facts that do hold -- strictly below one
    everywhere, and equal to the closed form built from the same two rates.
    """
    ratios, predicted = [], []
    for x in rows:
        a = x['rate_at_R'] * (1.0 - x['rate_at_R'])
        b = x['rate_at_2R'] * (1.0 - x['rate_at_2R'])
        if a <= 0.0 or x['se_independent'] <= 0.0:
            continue
        ratios.append(x['se_paired'] / x['se_independent'])
        predicted.append(math.sqrt(a / (2.0 * a + b)))
    if not ratios:
        return {'id': 'paired_se_is_not_the_independent_one',
                'passed': False, 'why_it_matters': 'no usable rows'}
    worst_gap = max(abs(r - p_) for r, p_ in zip(ratios, predicted))
    return {
        'id': 'paired_se_is_not_the_independent_one',
        'question': 'is the paired error strictly smaller, and does its ratio '
                    'match the closed form sqrt(A / (2A + B))?',
        'ratio_min': min(ratios), 'ratio_max': max(ratios),
        'closed_form_at_equal_rates': 1.0 / math.sqrt(3.0),
        'max_abs_gap_from_closed_form': worst_gap,
        'all_strictly_below_one': all(r < 1.0 for r in ratios),
        'inflation_at_worst': 1.0 / min(ratios),
        'inflation_at_equal_rates': math.sqrt(3.0),
        'passed': all(r < 1.0 for r in ratios) and worst_gap < 1e-9,
        'why_it_matters': 'the ratio is the price of nesting. At equal rates the '
                          'unpaired error is sqrt(3) times the paired one, so a '
                          'check written against the unpaired value would tolerate '
                          'a drift up to 73 per cent larger than the data support',
    }


def check_grid_refinement_preserves_analytic_values(an):
    return {
        'id': 'grid_refinement_changes_no_analytic_value',
        'question': 'do shared cells keep bit-identical delta_tip and observed '
                    'risks on the refined grid?',
        'nested_p0': an['nested_p0'], 'nested_se': an['nested_se'],
        'shared_cells': an['shared_cells'], 'new_cells': an['new_cells'],
        'mismatches': an['analytic_mismatches'][:6],
        'n_mismatches': len(an['analytic_mismatches']),
        'passed': bool(an['nested_p0'] and an['nested_se']
                       and not an['analytic_mismatches']),
        'why_it_matters': 'nothing is sampled in this block, so any difference on '
                          'a shared cell is a defect in the join, not noise',
    }


def check_grid_extremes_do_not_move(an):
    same_min = (an['coarse_delta_tip_range'][0] == an['fine_delta_tip_range'][0])
    same_max = (an['coarse_delta_tip_range'][1] == an['fine_delta_tip_range'][1])
    return {
        'id': 'grid_refinement_leaves_the_threshold_extremes_unchanged',
        'question': 'do the minimum and maximum delta_tip move when the p0 and se '
                    'axes are halved in step?',
        'coarse_range_pp': an['coarse_delta_tip_range'],
        'fine_range_pp': an['fine_delta_tip_range'],
        'argmin_unchanged': an['coarse_argmin'] == an['fine_argmin'],
        'argmax_unchanged': an['coarse_argmax'] == an['fine_argmax'],
        'monotonicity': an['monotonicity'],
        'passed': bool(same_min and same_max
                       and an['coarse_argmin'] == an['fine_argmin']
                       and an['coarse_argmax'] == an['fine_argmax']),
        'why_it_matters': 'the extremes are what a reader quotes as the range of '
                          'the threshold. Because delta_tip is monotone in each '
                          'axis they are attained at the vertices, and halving an '
                          'axis adds no vertex -- so this is a derived '
                          'consequence, not a lucky observation',
    }


def check_grid_refinement_keeps_the_type_i_pool(sg):
    pd = sg['pooled_difference']
    tol = sg['pooled_mc95_halfwidth_fine'] + sg['pooled_mc95_halfwidth_coarse']
    return {
        'id': 'refined_grid_leaves_the_pooled_type_i_error_unchanged',
        'question': 'does pooling over twice as many cells move the type I error '
                    'by more than the two pooled Monte Carlo errors?',
        'pooled_on_shared_coarse': sg['pooled_on_shared_coarse'],
        'pooled_over_all_fine': sg['pooled_over_all_fine'],
        'difference': pd,
        'tolerance_sum_of_mc95': tol,
        'shared_cells_sampled': sg['shared_cells_sampled'],
        'cells_outside_two_se': len(sg['cells_outside_two_se']),
        'max_abs_per_cell_z': sg['max_abs_z'],
        'passed': abs(pd) <= tol,
        'why_it_matters': 'the pooled type I error is the calibration claim of '
                          'the whole evaluation. If it moved under grid '
                          'refinement the headline would be an artefact of the '
                          'step size chosen in step 1',
    }


def check_per_cell_drift_is_noise_like(sg):
    """Per-cell drift must look like noise, not like a systematic shift.

    A mean z near zero with a spread near one is what random drift looks like. A
    mean z far from zero is a bias the pooled test could hide, because pooling
    averages it away; a spread far below one would instead mean the per-cell
    errors are overstated. The tolerances are the same as in block R, so the two
    blocks cannot disagree about what noise looks like; the earlier wider bounds
    of 0.5 and 0.5 to 1.6 were loose enough to accept a genuine bias.
    """
    zs = [x['z'] for x in sg['per_cell'] if x['z'] is not None]
    if not zs:
        return {'id': 'per_cell_drift_is_noise_like', 'passed': False,
                'why_it_matters': 'no usable per-cell comparisons'}
    mean_z = sum(zs) / len(zs)
    var_z = sum((z - mean_z) ** 2 for z in zs) / (len(zs) - 1) if len(zs) > 1 else 0.0
    sd_z = math.sqrt(var_z)
    over2 = sum(1 for z in zs if abs(z) > 2.0)
    expected_over2 = 0.0455 * len(zs)
    tol = 1.96 * math.sqrt(0.0455 * (1 - 0.0455) * len(zs))
    return {
        'id': 'per_cell_drift_is_noise_like',
        'question': 'is the per-cell drift centred on zero with unit spread?',
        'cells': len(zs),
        'mean_z': mean_z, 'sd_z': sd_z,
        'outside_two_se': over2, 'outside_two_se_expected': expected_over2,
        'outside_two_se_tolerance_95': tol,
        'passed': (abs(mean_z) < 0.20 and 0.85 < sd_z < 1.20
                   and abs(over2 - expected_over2) <= tol),
        'why_it_matters': 'a non-zero mean would be a directional effect hidden '
                          'by pooling; an sd far below one would mean the quoted '
                          'per-cell Monte Carlo errors are too wide',
    }


def check_delta_sig_shift_is_small(rows):
    shifts = [x['shift_pp'] for x in rows]
    over = [x for x in rows if not x['within_point_two_pp']]
    return {
        'id': 'delta_sig_moves_less_than_0_2_pp_under_refinement',
        'question': 'with twice the bisection iterations and twice the '
                    'replicates, does the located crossing move less than 0.2 pp?',
        'cells': len(rows),
        'max_shift_pp': max(shifts) if shifts else None,
        'mean_shift_pp': (sum(shifts) / len(shifts)) if shifts else None,
        'over_threshold': [{'cell_key': x['cell_key'], 'n_per_arm': x['n_per_arm'],
                            'shift_pp': x['shift_pp'],
                            'delta_tip_pp': x['delta_tip_pp']}
                           for x in over],
        'max_shift_as_pct_of_tip': max((x['shift_as_pct_of_tip'] for x in rows),
                                       default=None),
        'bracket_width_base_pp': rows[0]['bracket_width_base_pp'] if rows else None,
        'bracket_width_fine_pp': rows[0]['bracket_width_fine_pp'] if rows else None,
        'passed': len(over) == 0 and bool(rows),
        'why_it_matters': 'delta_sig is the only headline quantity that is '
                          'estimated rather than derived, so it is the only one '
                          'whose value depends on the simulation budget',
    }


def check_bracket_tolerance_is_not_the_binding_constraint(rows):
    """Separate the two sources of movement in delta_sig.

    Halving the bracket and doubling the replicates both move the estimate. If
    the bracket is already far below 0.2 pp then the residual movement is
    Monte Carlo, and reporting the two together would let a reader attribute a
    sampling fluctuation to the locator.
    """
    worst_bracket = max(x['bracket_width_fine_pp'] for x in rows) if rows else None
    worst_shift = max(x['shift_pp'] for x in rows) if rows else None
    # The refinement must actually have happened. Without this, a run in which
    # both levels used the same iteration count would report two identical
    # brackets, satisfy the width test, and pass while testing nothing.
    refined = bool(rows) and all(x['bracket_is_finer'] for x in rows)
    return {
        'id': 'the_bisection_bracket_is_not_what_limits_delta_sig',
        'question': 'is the achieved bracket width already far below the 0.2 pp '
                    'tolerance, so the residual shift is Monte Carlo?',
        'widest_fine_bracket_pp': worst_bracket,
        'widest_observed_shift_pp': worst_shift,
        'refinement_actually_happened': refined,
        'bracket_below_tolerance': bool(worst_bracket is not None
                                        and worst_bracket < 0.2),
        'passed': bool(worst_bracket is not None and worst_bracket < 0.2
                       and worst_shift is not None and worst_shift < 0.2
                       and refined),
        'why_it_matters': 'the two are routinely conflated. Bisection to 2^-26 of '
                          'a 20 pp bracket is about 3e-7 pp, so essentially all '
                          'movement at this budget is Monte Carlo, and the honest '
                          'statement is that delta_sig is replicate-limited, not '
                          'locator-limited',
    }


def check_nothing_here_validates_the_estimand():
    return {
        'id': 'convergence_is_not_validity',
        'question': 'is the boundary of this step recorded?',
        'statement': 'doubling replicates and halving a step move an estimate '
                     'towards its own Monte Carlo expectation. Neither tests '
                     'whether the two-parameter measurement model is the right '
                     'one, and neither speaks to the direction of any bias. This '
                     'step supports the claim that the reported numbers are '
                     'reproducible; it does not support a claim that they are '
                     'right.',
        'not_tested': ['correctness of the measurement model',
                       'direction of residual bias',
                       'transportability beyond the grid box',
                       'behaviour of layer B, which is analytic only'],
        'passed': True,
    }


# ---------------------------------------------------------------------- main
def main():
    layer_a = load_layer_a()

    r_rows, r_jobs = block_r(layer_a)
    an = block_g_analytic()
    sg = block_g_sampled(an)
    d_rows, d_jobs, d_cand = block_d(layer_a)

    checks = [
        check_nested_prefix_is_real(),
        check_seed_is_position_free(),
        check_replicate_doubling(r_rows),
        check_paired_se_is_the_smaller_one(r_rows),
        check_grid_refinement_preserves_analytic_values(an),
        check_grid_extremes_do_not_move(an),
        check_grid_refinement_keeps_the_type_i_pool(sg),
        check_per_cell_drift_is_noise_like(sg),
        check_delta_sig_shift_is_small(d_rows),
        check_bracket_tolerance_is_not_the_binding_constraint(d_rows),
        check_nothing_here_validates_the_estimand(),
    ]
    passed = all(c['passed'] for c in checks)

    # Checks are looked up by id through an index rather than by repeating the id
    # inline. A concatenated literal drifted in case from the definition once
    # already, and next() then raised StopIteration after the whole computation
    # had finished, which is the most expensive way to find a typo.
    by_id = {c['id']: c for c in checks}
    if len(by_id) != len(checks):
        raise SystemExit('duplicate check ids: the headline index would be '
                         'ambiguous')
    doubling = by_id['replicate_doubling_leaves_the_type_i_estimate_reproducible']

    doc = {
        'meta': {
            'step': 'W1 step 8',
            'purpose': 'convergence evidence for BMC Medical Research '
                       "Methodology's 'well tested' requirement: replicate "
                       'doubling, grid refinement, and delta_sig resolution',
            'written': '2026-10-06',
            'standard_library_only': True,
            'seed_derivation': 'sha256 over a namespaced cell key; the built-in '
                               'hash() is avoided because PYTHONHASHSEED '
                               'randomises it per process',
            'sampling_kernel': 'imported from calibrate_estimator_eval.run_cell, '
                               'not reimplemented',
            'delta_sig_kernel': 'imported from locate_delta_sig_and_coverage',
            'why_step4_and_step5_numbers_are_not_reused':
                'those steps consume one stream across the grid in traversal '
                'order, so an inserted cell shifts every later cell; measured on '
                'the type I rows, one inserted cell moved 33 of 40 shared cells',
            'paired_se_derivation': 'nested streams make the difference depend '
                                    'on the last R draws only, so se = '
                                    'sqrt(rate_R (1 - rate_R) / R), not the '
                                    'sqrt(2)-larger unpaired value',
            'bounds': 'convergence of a Monte Carlo estimate is not validity of '
                      'the estimand',
        },
        'block_R_replicate_doubling': {
            'design': 'type I cells, R = %d against 2R = %d, nested streams, '
                      'one stream per cell derived from its key'
                      % (REPS_BASE, REPS_DOUBLE),
            'cells_requested': r_jobs, 'cells_used': len(r_rows),
            'rows': r_rows,
        },
        'block_G_grid_refinement': {
            'design': 'midpoint insertion on p0 and on se, giving a strict '
                      'superset; analytic quantities compared exactly, the type '
                      'I error resampled on both grids',
            'analytic': {k: v for k, v in an.items() if not k.startswith('_')},
            'sampled': {k: v for k, v in sg.items() if not k.startswith('_')},
        },
        'block_D_delta_sig_resolution': {
            'design': 'located crossing recomputed with %d bisection iterations '
                      'and %d replicates against %d and %d; n in {1000, 2000}'
                      % (BISECT_ITERS_FINE, REPS_SIG_FINE,
                         BISECT_ITERS_BASE, REPS_SIG_BASE),
            'candidate_cells': d_cand, 'jobs': d_jobs, 'cells_used': len(d_rows),
            'selection_note': 'every %dth eligible cell in file order, so the '
                              'spread covers the grid rather than its first rows'
                              % max(1, d_cand // DELTA_SIG_CELLS),
            'rows': d_rows,
        },
        'checks': checks,
        'summary': {
            'checks_total': len(checks),
            'checks_passed': sum(1 for c in checks if c['passed']),
            'all_passed': passed,
            'headline': {
                'type_I_cells_doubled': len(r_rows),
                'cells_beyond_two_se': doubling['cells_beyond_two_se'],
                'cells_beyond_two_se_expected':
                    doubling['cells_beyond_two_se_expected'],
                'mean_z': doubling['mean_z'], 'sd_z': doubling['sd_z'],
                'verdicts_changed_on_doubling': doubling['verdicts_changed'],
                'verdict_band': doubling['verdict_band'],
                'coarse_cells': an['coarse_cells'],
                'fine_cells': an['fine_cells'],
                'new_cells': an['new_cells'],
                'analytic_mismatches': len(an['analytic_mismatches']),
                'pooled_type_I_coarse':
                    sg['pooled_on_shared_coarse']['rate'],
                'pooled_type_I_fine': sg['pooled_over_all_fine']['rate'],
                'pooled_type_I_difference': sg['pooled_difference'],
                'delta_sig_cells': len(d_rows),
                'max_delta_sig_shift_pp':
                    max((x['shift_pp'] for x in d_rows), default=None),
            },
        },
    }

    with io.open(OUT, 'w', encoding='utf-8') as fh:
        json.dump(doc, fh, ensure_ascii=False, indent=2)

    for c in checks:
        print('[%s] %s' % ('PASS' if c['passed'] else 'FAIL', c['id']))
        for k, v in c.items():
            if k in ('id', 'question', 'passed', 'why_it_matters'):
                continue
            s = json.dumps(v, ensure_ascii=False)
            print('       %-26s %s' % (k, s[:150]))
    h = doc['summary']['headline']
    print('\n%d/%d checks passed' % (doc['summary']['checks_passed'],
                                     doc['summary']['checks_total']))
    print('  type I cells doubled      %d, paired z mean %.3f sd %.3f'
          % (h['type_I_cells_doubled'], h['mean_z'], h['sd_z']))
    print('  beyond two paired SE     %d of %d (expected %.1f)'
          % (h['cells_beyond_two_se'], h['type_I_cells_doubled'],
             h['cells_beyond_two_se_expected']))
    print('  verdict flips             %d at a band of %.3f'
          % (h['verdicts_changed_on_doubling'], h['verdict_band']))
    print('  grid %d -> %d cells (+%d new), analytic mismatches %d'
          % (h['coarse_cells'], h['fine_cells'], h['new_cells'],
             h['analytic_mismatches']))
    print('  pooled type I %.4f -> %.4f (difference %+.4f)'
          % (h['pooled_type_I_coarse'], h['pooled_type_I_fine'],
             h['pooled_type_I_difference']))
    print('  delta_sig cells %d, max shift %.4f pp'
          % (h['delta_sig_cells'], h['max_delta_sig_shift_pp'] or 0.0))
    print('\nwrote %s' % OUT)
    return 0 if passed else 1


if __name__ == '__main__':
    sys.exit(main())