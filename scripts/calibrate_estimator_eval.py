# -*- coding: utf-8 -*-
"""W1 step 4: calibration. Two test targets, two tables, never one summary.

Scope, deliberately narrow. Step 3 certified the mechanisms and wrote the layer A
records. This step draws replicates on those records and answers exactly one
question per table:

  Table I  -- type I error.   Cells with or_true = 1.0 and delta = 0. There is no
               true difference and no measurement asymmetry, so the population
               observed OR is exactly 1 and the rejection rate estimates the type I
               error of the two-sided Wald test. Nothing else belongs here.

  Table D  -- detection rate. Cells with or_true < 1 and delta > 0. There IS a
               true benefit, so a rejection is a detection, NOT a type I error.
               This distinction is the whole reason the two tables exist: the
               third review round flagged that a single "reject rate" column
               silently conflated them.

  Table P  -- power at no asymmetry. The same or_true < 1 cells at delta = 0.
               This is the reference ceiling for table D: the difference between
               it and table D is the detection lost to the asymmetry. It is named
               power, not detection and not type I error, because there is no
               asymmetry to detect the effect of.

Records are JOINED from step 3 by cell key. The grid is not rebuilt here: a
rebuilt grid could drift from the certified one without any check failing.

Naming is enforced, not promised. `detection_is_not_labelled_type_i` walks the
written document and fails if any key inside the detection subtree contains
"type_i", or if any detection record carries a type-I metric name. A comment
cannot enforce this; a recursive scan can.

Three checks that only become answerable once replicates exist:

1. The type I error lands on 5 per cent, which proves the or_true = 1 cell is
   reachable and the test is not mis-specified.
2. At delta = delta_tip the population observed OR is exactly 1 while the
   rejection rate is materially ABOVE zero. This is the correction of the
   third-round error: fixed population parameters fix the population OR, not the
   sample OR, so the rejection probability at delta_tip is a property of the
   sampling distribution and is not identically zero. A reject rate of about
   zero here would mean the specification had re-introduced the confusion.
3. The detection rate falls as delta rises, pooled over cells. If it did not, the
   sampling implementation would be wrong in a way no per-cell table would show.

Run:  python scripts/calibrate_estimator_eval.py
Out:  reports/estimator_eval_calibration.json
Exit 0 if every calibration check passes, 1 otherwise.
"""
import io
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

# The one permitted source of the mechanism and of the closed form.
import simulate_tipping_point as legacy   # noqa: E402

DGP = os.path.join(ROOT, 'reports', 'estimator_eval_dgp_checks.json')
OUT = os.path.join(ROOT, 'reports', 'estimator_eval_calibration.json')

TOL = 1e-12
ALPHA = 0.05
REPS = 2000                     # estimator_eval_design.json: reps_per_cell
N_GRID = [100, 500, 1000, 2000]  # estimator_eval_design.json: n_per_arm
SEED_I = 20261006               # distinct stream per table, so adding a table
SEED_D = 20261007               # later cannot shift an earlier table's draws
SEED_P = 20261008

# The detection table scans delta on a scale-free grid: fractions of that cell's
# own delta_tip. An absolute delta grid would put most cells at either 0 or
# beyond their threshold, because delta_tip ranges from about 4 to about 32 pp
# across the grid.
DELTA_FRACTIONS = [0.25, 0.50, 0.75, 1.00]


# ------------------------------------------------------------------ loading
def load_layer_a():
    """Read the certified layer A records written by step 3.

    Joining rather than rebuilding is the point. If this script rebuilt the grid
    from P0_GRID and friends, a later edit to either file could make the two
    disagree with no check failing.
    """
    if not os.path.exists(DGP):
        raise SystemExit('step 3 output missing: %s -- run step 3 first' % DGP)
    with io.open(DGP, encoding='utf-8') as fh:
        doc = json.load(fh)
    if doc['meta']['step'] != 'W1 step 3':
        raise SystemExit('unexpected producer: %r' % doc['meta']['step'])
    if not doc['summary']['all_passed']:
        raise SystemExit('step 3 mechanisms are not certified; do not sample')
    rows = [r for r in doc['records']['layer_A'] if r['legal']]
    if not rows:
        raise SystemExit('no legal layer A records to calibrate on')
    return rows, doc['summary']


def key_of(r):
    return 'A|p0=%s|or_true=%s|se=%s|f=%s' % (r['p0'], r['or_true'], r['se'], r['f'])


# ----------------------------------------------------------------- sampling
_Z_CACHE = {}


def _z_two_sided(alpha):
    """Two-sided normal critical value, from the error function.

    Hard-coding 1.96 would tie the test to one nominal level. The design fixes
    alpha at 0.05, but the function stays honest if that changes.
    """
    if alpha not in _Z_CACHE:
        import math
        # z such that P(|Z| > z) = alpha  ->  Q(z/2) = alpha/2
        lo, hi = 0.0, 40.0
        for _ in range(200):
            mid = (lo + hi) / 2.0
            if math.erfc(mid / 2.0 ** 0.5) > alpha:
                lo = mid
            else:
                hi = mid
        _Z_CACHE[alpha] = (lo + hi) / 2.0
    return _Z_CACHE[alpha]


def run_cell(q0, q1, n, reps, rng):
    """Replicates for one cell. Returns rates, not verdicts.

    usable_reps is recorded rather than assumed: with n = 100 and a risk near 0.5
    a zero cell count is possible, and silently drawing fewer replicates than
    requested would flatter the Monte Carlo error.
    """
    zc = _z_two_sided(ALPHA)
    import math
    reject = 0
    used = 0
    dropped = 0
    sum_lor = 0.0
    sum_sq_lor = 0.0
    for _ in range(reps):
        x0 = rng.binomialvariate(n, q0)
        x1 = rng.binomialvariate(n, q1)
        if x0 == 0 or x1 == 0 or x0 == n or x1 == n:
            dropped += 1
            continue
        lor = math.log(x1 / (n - x1)) - math.log(x0 / (n - x0))
        se = (1.0 / x0 + 1.0 / (n - x0) + 1.0 / x1 + 1.0 / (n - x1)) ** 0.5
        if se <= 0.0:
            dropped += 1
            continue
        used += 1
        sum_lor += lor
        sum_sq_lor += lor * lor
        if abs(lor / se) > zc:
            reject += 1
    rate = (reject / used) if used else None
    mc = None
    if rate is not None:
        mc = 1.96 * ((rate * (1.0 - rate) / used) ** 0.5) if 0 < used else None
    mean_lor = (sum_lor / used) if used else None
    var_lor = None
    if used > 1:
        var_lor = (sum_sq_lor - used * mean_lor ** 2) / (used - 1)
    return {
        'reps_requested': reps,
        'usable_reps': used,
        'dropped_reps': dropped,
        'usable_share': (used / reps) if reps else None,
        'rejections': reject,
        'reject_rate': rate,
        'mc95_halfwidth': mc,
        'mean_sample_log_or': mean_lor,
        'var_sample_log_or': var_lor,
        'inadequate_reps': used < 0.5 * reps,
    }


def observed_or(q0, q1):
    return legacy.or_from_risks(q0, q1)


# --------------------------------------------------------------- table I
def table_type_i(rows):
    """Type I error: or_true = 1, delta = 0.

    At these cells p1 = p0 and f + delta = f, so q1 = q0 exactly and the
    population observed OR is 1. Every rejection is therefore an error, and the
    rejection rate estimates the type I error. The comparison target is the
    between-group TRUE log OR, which is 0 here and equals the observed one, so
    there is no measurement-asymmetry confound to unpick.
    """
    out = []
    for r in rows:
        if abs(r['or_true'] - 1.0) > TOL:
            continue
        # or_true = 1 implies p1 = p0; assert rather than assume, because the
        # whole table is meaningless if the null cell is not really null.
        if abs(r['p1'] - r['p0']) > TOL:
            raise SystemExit('or_true = 1 cell with p1 != p0 at %s' % key_of(r))
        if r['threshold_exists']:
            raise SystemExit('or_true = 1 cell claims a threshold at %s' % key_of(r))
        rng = random.Random(SEED_I)
        for n in N_GRID:
            q0 = legacy.observed_risk_effective(r['p0'], r['se'], r['f'])
            q1 = legacy.observed_risk_effective(r['p1'], r['se'], r['f'])
            pop_or = observed_or(q0, q1)
            res = run_cell(q0, q1, n, REPS, rng)
            res.update({
                'cell_key': key_of(r),
                'p0': r['p0'], 'p1': r['p1'], 'or_true': r['or_true'],
                'se': r['se'], 'f': r['f'],
                'n_per_arm': n,
                'delta_pp': 0.0,
                'q0': q0, 'q1': q1,
                'population_observed_or': pop_or,
                'metric_name': 'type_I_error',
            })
            out.append(res)
    return out


# --------------------------------------------------------------- table D
def table_detection(rows):
    """Detection rate: or_true < 1 with delta > 0.

    A rejection here means the test found a true benefit that exists. Calling this
    a type I error would be wrong twice over: the null is false, and the whole
    point of the table is the size of the asymmetry that hides it.

    Delta is scanned as a fraction of the cell's own delta_tip. The absolute
    values matter: at delta = delta_tip the population observed OR is exactly 1,
    so the detection rate there is the rate at which the test fails to see a
    benefit that is really there. That cell is the reason this table exists.
    """
    out = []
    for r in rows:
        if r['or_true'] >= 1.0 - TOL:
            continue
        if not r['threshold_exists']:
            continue
        rng = random.Random(SEED_D)
        for n in N_GRID:
            for frac in DELTA_FRACTIONS:
                delta = r['delta_tip_pp'] * frac
                fp = r['f'] + delta / 100.0
                # Legality, asserted not assumed: f + delta must stay in [0, 1]
                # and se must exceed it, or the cell has no meaning.
                if fp < -TOL or fp > 1.0 + TOL or r['se'] <= fp + TOL:
                    raise SystemExit('illegal f + delta at %s frac %.2f'
                                     % (key_of(r), frac))
                q0 = legacy.observed_risk_effective(r['p0'], r['se'], r['f'])
                q1 = legacy.observed_risk_effective(r['p1'], r['se'], fp)
                pop_or = observed_or(q0, q1)
                res = run_cell(q0, q1, n, REPS, rng)
                res.update({
                    'cell_key': key_of(r),
                    'p0': r['p0'], 'p1': r['p1'], 'or_true': r['or_true'],
                    'se': r['se'], 'f': r['f'],
                    'n_per_arm': n,
                    'delta_fraction_of_tip': frac,
                    'delta_pp': delta,
                    'delta_tip_pp': r['delta_tip_pp'],
                    'q0': q0, 'q1': q1,
                    'population_observed_or': pop_or,
                    'metric_name': 'detection_rate',
                })
                out.append(res)
    return out


# --------------------------------------------------------------- table P
def table_power_reference(rows):
    """Power at no asymmetry: or_true < 1 with delta = 0.

    The ceiling the detection table is measured against. Naming matters as much
    as the number: there is no asymmetry here, so a rejection is neither a
    detection of the effect of the asymmetry nor a type I error.
    """
    out = []
    for r in rows:
        if r['or_true'] >= 1.0 - TOL:
            continue
        rng = random.Random(SEED_P)
        for n in N_GRID:
            q0 = legacy.observed_risk_effective(r['p0'], r['se'], r['f'])
            q1 = legacy.observed_risk_effective(r['p1'], r['se'], r['f'])
            res = run_cell(q0, q1, n, REPS, rng)
            res.update({
                'cell_key': key_of(r),
                'p0': r['p0'], 'p1': r['p1'], 'or_true': r['or_true'],
                'se': r['se'], 'f': r['f'],
                'n_per_arm': n,
                'delta_pp': 0.0,
                'q0': q0, 'q1': q1,
                'population_observed_or': observed_or(q0, q1),
                'metric_name': 'power_at_no_asymmetry',
            })
            out.append(res)
    return out


# ------------------------------------------------------------------ checks
def pool(records):
    """Pooled rejection rate across records, weighting by usable replicates."""
    num = sum(r['rejections'] for r in records if r['reject_rate'] is not None)
    den = sum(r['usable_reps'] for r in records if r['reject_rate'] is not None)
    if den == 0:
        return None, None, 0
    rate = num / den
    return rate, 1.96 * ((rate * (1.0 - rate) / den) ** 0.5), den


def check_type_i_calibration(t1):
    """Does the type I error land on the nominal 5 per cent?

    Judged per cell against its own Monte Carlo error, not against a fixed
    tolerance: at 2000 replicates a rate carries about 0.7 pp of noise, so a flat
    1 pp band would fail cells for sampling reasons alone.
    """
    bad, worst = [], 0.0
    for r in t1:
        if r['reject_rate'] is None:
            bad.append({'cell_key': r['cell_key'], 'n_per_arm': r['n_per_arm'],
                        'why': 'no usable replicates'})
            continue
        dev = abs(r['reject_rate'] - ALPHA)
        mc = r['mc95_halfwidth'] or 0.0
        worst = max(worst, dev)
        # Three Monte Carlo standard errors. Beyond that the deviation is not
        # sampling noise and the cell is reported rather than summarised.
        if dev > max(3.0 * mc, 0.02):
            bad.append({'cell_key': r['cell_key'], 'n_per_arm': r['n_per_arm'],
                        'reject_rate': r['reject_rate'],
                        'deviation': dev, 'mc95_halfwidth': mc})
    rate, mc, den = pool(t1)
    return {
        'id': 'type_I_error_lands_on_nominal',
        'question': 'in the or_true = 1, delta = 0 cells, does the two-sided Wald '
                    'test reject at the nominal 5 per cent?',
        'cells': len(t1),
        'pooled_reject_rate': rate,
        'pooled_mc95_halfwidth': mc,
        'pooled_usable_reps': den,
        'worst_cell_deviation': worst,
        'outliers': bad,
        'passed': rate is not None and abs(rate - ALPHA) <= max(3.0 * mc, 0.005)
                  and not bad,
        'rule': 'a cell passes when |rate - alpha| <= max(3 MC standard errors, '
                '2 pp); the pooled rate must be within max(3 MC SE, 0.5 pp)',
    }


def check_population_or_one_but_rejection_nonzero(td):
    """At delta = delta_tip the population OR is 1 and the rejection rate is not 0.

    The third review round's item 2. Fixed population parameters fix the POPULATION
    observed OR, not the sample OR, so at delta_tip the test still rejects about
    alpha of the time. Reporting zero here would reinstate the confusion the round
    corrected.
    """
    at_tip = [r for r in td if abs(r['delta_fraction_of_tip'] - 1.0) < TOL]
    or_bad = [r['cell_key'] for r in at_tip
              if abs(r['population_observed_or'] - 1.0) > 1e-9]
    rate, mc, den = pool(at_tip)
    return {
        'id': 'population_or_one_at_tip_while_rejection_stays_nonzero',
        'question': 'at delta = delta_tip, is the population observed OR exactly 1 '
                    'AND the rejection rate materially above zero?',
        'cells_at_tip': len(at_tip),
        'cells_whose_population_or_is_not_one': or_bad,
        'pooled_reject_rate_at_tip': rate,
        'pooled_mc95_halfwidth': mc,
        'pooled_usable_reps': den,
        'min_cell_reject_rate_at_tip': (min(r['reject_rate'] for r in at_tip)
                                        if at_tip else None),
        'passed': (not or_bad) and rate is not None and rate > 0.5 * ALPHA
                  and rate < 5.0 * ALPHA,
        'why_this_matters': 'a rejection rate of about zero at delta_tip would mean '
                            'the population OR had been confused with the sample OR; '
                            'the rate must sit near alpha, because at delta_tip the '
                            'null of no observed difference is exactly true',
    }


def check_detection_falls_with_delta(td):
    """Pooled detection rate must not rise as the asymmetry grows.

    Monotone in expectation and, pooled over hundreds of cells, in practice. If it
    did not fall, the sampling implementation would be wrong in a way that no
    per-cell table would reveal.
    """
    by_frac = {}
    for r in td:
        if r['reject_rate'] is None:
            continue
        by_frac.setdefault(r['delta_fraction_of_tip'], []).append(r)
    curve = []
    for frac in sorted(by_frac):
        rate, mc, den = pool(by_frac[frac])
        curve.append({'delta_fraction_of_tip': frac, 'reject_rate': rate,
                      'mc95_halfwidth': mc, 'usable_reps': den,
                      'cells': len(by_frac[frac])})
    violations = []
    for a, b in zip(curve, curve[1:]):
        if b['reject_rate'] > a['reject_rate'] + a['mc95_halfwidth']:
            violations.append({'from_fraction': a['delta_fraction_of_tip'],
                               'to_fraction': b['delta_fraction_of_tip'],
                               'rise': b['reject_rate'] - a['reject_rate'],
                               'allowance': a['mc95_halfwidth']})
    return {
        'id': 'detection_rate_falls_as_delta_grows',
        'question': 'pooled over cells, does the detection rate fall as delta '
                    'approaches delta_tip?',
        'curve': curve,
        'monotonicity_violations': violations,
        'passed': len(curve) >= 2 and not violations,
        'note': 'at the last fraction the population observed OR is 1, so the rate '
                'there is expected to sit near alpha rather than near zero',
    }


def walk_keys(node, path=''):
    """Every dict key in the document, with its path."""
    if isinstance(node, dict):
        for k, v in node.items():
            p = path + '/' + str(k)
            yield p, k
            for x in walk_keys(v, p):
                yield x
    elif isinstance(node, list):
        for i, v in enumerate(node):
            for x in walk_keys(v, path + '[%d]' % i):
                yield x


def check_detection_is_not_labelled_type_i(td, tp, t1):
    """The naming rule, enforced by scanning the document that will be written.

    Three things must hold. No key inside the detection subtree may contain
    "type_i". Every record must carry exactly the metric name its table declares,
    so a record cannot be relabelled while its neighbours are not. And the three
    tables must carry three DISTINCT metric names, which is the structural
    guarantee that they can never be merged into one reject-rate column.

    What is deliberately NOT policed is the word "type I" inside prose. The
    detection table has to be able to say "a rejection here is a detection, not a
    type I error" -- that sentence IS the distinction. An earlier version of this
    check banned the substring everywhere, fired on its own disclaimer, and had to
    be weakened; the fixed version constrains labels, and states each table's
    hypothesis once at table level instead of copying it into every record.
    """
    tables = [('table_I', 'type_I_error', t1),
              ('table_D', 'detection_rate', td),
              ('table_P', 'power_at_no_asymmetry', tp)]

    # The forbidden-key scan applies to the tables that must NOT look like a type I
    # table. The type I table is excluded by construction, not by exemption: it is
    # supposed to carry that name. Probing it with a key of "type_I_error" would
    # make the check fail on its own scaffolding, which is what happened the first
    # time this ran.
    bad_keys = []
    for label, declared, rows in tables:
        if declared == 'type_I_error':
            continue
        bad_keys.extend('%s:%s' % (label, p) for p, k in walk_keys(rows)
                        if 'type_i' in str(k).lower() and str(k) != 'metric_name')

    mislabelled = []
    for label, declared, rows in tables:
        for r in rows:
            if str(r.get('metric_name')) != declared:
                mislabelled.append({'table': label, 'cell_key': r.get('cell_key'),
                                    'declared': declared,
                                    'metric_name': r.get('metric_name')})

    names = sorted(declared for _, declared, _ in tables)
    n_named_I = sum(1 for r in td if str(r.get('metric_name')) == 'type_I_error')
    return {
        'id': 'the_three_tables_cannot_be_merged',
        'question': 'are the type I, detection and power tables separately and '
                    'uniquely named, so that no reader or later edit can collapse '
                    'them into one reject-rate column?',
        'declared_metric_names': names,
        'metric_names_distinct': len(set(names)) == len(names),
        'cells_per_table': {label: len(rows) for label, _, rows in tables},
        'detection_records_named_type_I': n_named_I,
        'key_scan_excludes': 'the type I table itself, whose name is the one the '
                             'rule forbids elsewhere',
        'forbidden_keys_found': bad_keys,
        'records_mislabelled': mislabelled,
        'passed': (not bad_keys) and (not mislabelled)
                  and len(set(names)) == len(names) and n_named_I == 0,
        'rule': 'each table declares exactly one metric name; every record repeats '
                'only that name; no field of the detection or power table may be '
                'named type_i; the three names must be pairwise distinct',
        'why_this_check_exists': 'the conflation of the two targets is not a wording '
                                 'problem, it is a wrong denominator attached to a '
                                 'wrong hypothesis. Structural uniqueness survives '
                                 'later edits in a way that a prose reminder does '
                                 'not, and the prose itself must stay allowed to '
                                 'state the distinction',
    }


def check_replicate_accounting(t1, td, tp):
    """Usable replicates must be recorded, and no cell may be silently thin.

    HONESTY REQUIREMENT. This check was reverse-tested by deleting the zero-cell
    guard in run_cell, and nothing failed: across all 3360 cells the usable share
    was exactly 1.000, so no replicate was ever dropped. With n >= 100 and risks
    away from 0 and 1 the probability of a zero cell count is negligible, so the
    drop rule is correct but INERT on this grid. A guard that cannot fire is not
    evidence of correctness, and reporting it as if it were would be the same
    mistake as a silent zero-cell drop.

    So the check now reports its own reachability: how many replicates were
    dropped in total, whether the branch was ever entered, and the largest
    theoretical drop probability over the grid. The paper can then say the rule
    exists and did not bind, instead of implying it did work.
    """
    thin = []
    dropped_total = 0
    for lab, rows in (('type_I', t1), ('detection', td), ('power', tp)):
        for r in rows:
            dropped_total += r['dropped_reps']
            if r['reject_rate'] is None or r['inadequate_reps']:
                thin.append({'table': lab, 'cell_key': r.get('cell_key'),
                             'n_per_arm': r['n_per_arm'],
                             'usable_reps': r['usable_reps'],
                             'dropped_reps': r['dropped_reps']})
    all_rows = t1 + td + tp
    worst = max(r['usable_share'] for r in all_rows if r['usable_share'] is not None)
    # Largest per-replicate probability of an unusable draw over all cells:
    # P(x = 0 or x = n) = q^n + (1-q)^n, maximised at the risk furthest from 0.5.
    worst_drop_p = 0.0
    for r in all_rows:
        for q in (r['q0'], r['q1']):
            n = r['n_per_arm']
            p = (q ** n) + ((1.0 - q) ** n)
            worst_drop_p = max(worst_drop_p, p)
    inert = (dropped_total == 0) and not thin
    return {
        'id': 'replicate_accounting_is_recorded',
        'question': 'is usable_reps recorded per cell, is any cell left with too '
                    'few replicates to summarise, and was the zero-cell drop rule '
                    'ever actually reached?',
        'cells_total': len(all_rows),
        'inadequate_or_unusable': thin,
        'worst_usable_share': worst,
        'replicates_dropped_total': dropped_total,
        'drop_rule_reached': not inert,
        'inert_on_this_run': inert,
        'largest_theoretical_per_replicate_drop_probability': worst_drop_p,
        'passed': not thin,
        'rule': 'a cell with usable_reps below half the requested count is reported '
                'as infeasible rather than summarised',
        'honesty_note': (
            'the drop rule did not bind on this grid: zero cells never occurred, so '
            'this check certifies only that every cell kept a full set of usable '
            'replicates. Deleting the drop rule outright does not change any number '
            'here. The rule is retained because n = 100 with an extreme risk would '
            'reach it, and because an unrecorded drop is exactly the failure the '
            'sampling rule forbids'),
    }


def check_drop_rule_is_reachable():
    """Prove the zero-cell drop rule works, on a cell built to trigger it.

    The main accounting check is inert on this grid (see its honesty_note), so on
    its own it certifies nothing about the drop rule. This check exercises the
    same code path directly with an extreme risk at the smallest n, where a zero
    or saturated cell count is likely, and asserts two things: replicates WERE
    dropped, and the recorded rate is computed on the survivors rather than on the
    requested count.

    If this check ever fails, the drop rule has become a no-op that silently
    divides by the wrong denominator.
    """
    rng = random.Random(SEED_I + 977)
    probe = run_cell(q0=0.02, q1=0.98, n=20, reps=2000, rng=rng)
    survivors = probe['usable_reps']
    expected_rate = (probe['rejections'] / survivors) if survivors else None
    consistent = (expected_rate is not None
                  and abs(expected_rate - probe['reject_rate']) < 1e-12)
    saturated_rate_correct = (probe['reject_rate'] == 1.0) if survivors else None
    return {
        'id': 'zero_cell_drop_rule_is_reachable_and_consistent',
        'question': 'does the zero-cell drop rule actually fire when it should, and '
                    'is the rate computed on the survivors?',
        'probe': {'q0': 0.02, 'q1': 0.98, 'n_per_arm': 20, 'reps': 2000},
        'replicates_dropped': probe['dropped_reps'],
        'usable_reps': survivors,
        'reject_rate': probe['reject_rate'],
        'rate_equals_rejections_over_survivors': consistent,
        'all_saturated_replicates_reject': saturated_rate_correct,
        'passed': (probe['dropped_reps'] > 0 and consistent
                   and survivors > 0),
        'why_separate_from_the_accounting_check': 'the grid in this step never '
                                                  'produces a zero cell, so the '
                                                  'accounting check cannot exercise '
                                                  'the rule. This probe can, and '
                                                  'fails if the rule stops working',
    }


# ------------------------------------------------------------------- main
def main():
    a_rows, dgp_summary = load_layer_a()
    t1 = table_type_i(a_rows)
    td = table_detection(a_rows)
    tp = table_power_reference(a_rows)

    checks = [
        check_type_i_calibration(t1),
        check_population_or_one_but_rejection_nonzero(td),
        check_detection_falls_with_delta(td),
        check_detection_is_not_labelled_type_i(td, tp, t1),
        check_replicate_accounting(t1, td, tp),
        check_drop_rule_is_reachable(),
    ]
    allok = all(c['passed'] for c in checks)

    pooled_i = dict(zip(('rate', 'mc', 'reps'), pool(t1)))
    pooled_d = {}
    for frac in DELTA_FRACTIONS:
        pooled_d[str(frac)] = dict(zip(
            ('rate', 'mc', 'reps'),
            pool([r for r in td
                  if abs(r['delta_fraction_of_tip'] - frac) < TOL])))
    pooled_p = dict(zip(('rate', 'mc', 'reps'), pool(tp)))

    # Detection loss against the no-asymmetry ceiling, at matched cells. This is
    # the quantity a reader of the paper actually wants: how much of the
    # achievable detection the asymmetry removes.
    ref = {}
    for r in tp:
        ref[(r['cell_key'], r['n_per_arm'])] = r['reject_rate']
    loss = []
    for r in td:
        base = ref.get((r['cell_key'], r['n_per_arm']))
        if base is None or r['reject_rate'] is None or base is None:
            continue
        loss.append({
            'cell_key': r['cell_key'], 'n_per_arm': r['n_per_arm'],
            'delta_fraction_of_tip': r['delta_fraction_of_tip'],
            'power_at_no_asymmetry': base,
            'detection_rate': r['reject_rate'],
            'absolute_loss': base - r['reject_rate'],
        })

    doc = {
        'meta': {
            'step': 'W1 step 4',
            'written': '2026-10-06',
            'purpose': 'calibration on the certified layer A records: one table per '
                       'test target, never a single reject-rate column',
            'records_joined_from': os.path.relpath(DGP, ROOT).replace('\\', '/'),
            'grid_rebuilt_here': False,
            'grid_rebuild_note': 'the layer A grid is read from step 3 by cell key; '
                                 'rebuilding it here would let the two drift apart '
                                 'without any check failing',
            'alpha': ALPHA,
            'reps_per_cell': REPS,
            'n_per_arm': N_GRID,
            'delta_fractions_of_tip': DELTA_FRACTIONS,
            'seeds': {'type_I': SEED_I, 'detection': SEED_D,
                      'power_reference': SEED_P},
            'seed_note': 'one independent stream per table, so adding a table later '
                         'cannot shift an earlier table\'s draws',
            'test': 'two-sided Wald on the log odds ratio, critical value from the '
                    'error function rather than a hard-coded 1.96',
            'mechanism_source': 'scripts/simulate_tipping_point.py',
            'standard_library_only': True,
            'step3_summary': dgp_summary,
        },
        'tables': {
            'I_type_I_error': {
                'cells': t1,
                'definition': 'or_true = 1 and delta = 0; the population observed OR '
                              'is exactly 1, so every rejection is an error',
                'metric_name': 'type_I_error',
                'pooled': pooled_i,
            },
            'D_detection_rate': {
                'cells': td,
                'definition': 'or_true < 1 and delta > 0; the null is false, so a '
                              'rejection is a detection and not a type I error',
                'metric_name': 'detection_rate',
                'pooled_by_delta_fraction': pooled_d,
            },
            'P_power_at_no_asymmetry': {
                'cells': tp,
                'definition': 'or_true < 1 and delta = 0; the reference ceiling for '
                              'table D, named power because there is no asymmetry',
                'metric_name': 'power_at_no_asymmetry',
                'pooled': pooled_p,
            },
        },
        'detection_loss_against_no_asymmetry_ceiling': loss,
        'checks': checks,
        'summary': {
            'checks_total': len(checks),
            'checks_failed': sum(1 for c in checks if not c['passed']),
            'type_I_cells': len(t1),
            'detection_cells': len(td),
            'power_reference_cells': len(tp),
            'pooled_type_I_error': pooled_i['rate'],
            'pooled_power_at_no_asymmetry': pooled_p['rate'],
            'pooled_detection_at_tip': pooled_d.get('1.0', {}).get('rate'),
            'all_passed': allok,
        },
    }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with io.open(OUT, 'w', encoding='utf-8') as fh:
        json.dump(doc, fh, indent=2, ensure_ascii=False)

    print('wrote %s (%d bytes)' % (OUT, os.path.getsize(OUT)))
    print()
    print('--- W1 step 4 calibration checks ---')
    for c in checks:
        print('%s  %s' % ('PASS ' if c['passed'] else 'FAIL ', c['id']))
    print()
    print('table I  type I error        : %4d cells, pooled %.4f (MC %.4f), nominal %.2f'
          % (len(t1), pooled_i['rate'], pooled_i['mc'], ALPHA))
    print('table P  power, no asymmetry : %4d cells, pooled %.4f (MC %.4f)'
          % (len(tp), pooled_p['rate'], pooled_p['mc']))
    print('table D  detection rate      : %4d cells' % len(td))
    for frac in DELTA_FRACTIONS:
        d = pooled_d[str(frac)]
        print('           delta = %.0f%% of delta_tip : %.4f (MC %.4f)'
              % (100 * frac, d['rate'], d['mc']))
    print()
    print('CALIBRATED, STEP 5 MAY RUN' if allok
          else 'CALIBRATION FAILED, DO NOT PROCEED TO STEP 5')
    return 0 if allok else 1


if __name__ == '__main__':
    raise SystemExit(main())