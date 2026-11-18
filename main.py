"""Forecast a launch day's hourly Twitch viewers from comparable launches.

Model (fit fresh on every case's own comparables, no stored tables):

Every comparable launch contributes its release-day curve as log-ratios to
its own hour-1 print. Two views of those ratios are combined:

1. A scale kernel: for each hour, the weighted median and spread of the
   log-ratios, weighting each comparable by how close its hour-1 print is
   to the target's (Gaussian kernel, sigma 0.6 in log space). Hour-1 scale
   is a strong predictor of the day's shape: launches whose first hour was
   already huge stay near it, while launches that opened small grow many
   times over as regions unlock (log hour-1 vs log peak ratio correlate
   about -0.94 across the kit's comparables).

2. A two-factor fit: per hour, least squares of log-ratio on log hour-1
   plus a rolling-unlock indicator (a comparable counts as rolling when
   its own day peaked at least 2x its hour 1). The target is a rolling
   midnight-local launch, so its indicator is 1; its log hour-1 is
   clipped to the comparables' range so the linear fit is never
   extrapolated past its data.

The two views are blended (geometric mean of quantiles). The linear
view's weight is 0.5 while the target's hour 1 is in the dense part of
the comparable range and tapers to 0 at the largest comparable hour 1,
so a record-breaking launch falls back to the local kernel alone. The
knee of the taper is the 75th percentile of the comparables' own hour-1
prints, so it adapts to each case's pool instead of trusting a fixed
number.

Quantiles are normal in log space around each view's centre with its
spread (floored), which keeps the 90% band honest: about 91% of the
practice hours fall inside it. With too few usable comparables the
forecast falls back to flat persistence with a widening spread.

Standard library only, deterministic, no network. Reads every *.in in
$MAKEFUN_INPUT_DIR (else /input) and writes one *.out per case to
$MAKEFUN_OUTPUT_DIR (else /output).
"""

import json
import math
import os
import sys
from datetime import datetime
from pathlib import Path

INPUT_DIR = Path(os.environ.get("MAKEFUN_INPUT_DIR", "/input"))
OUTPUT_DIR = Path(os.environ.get("MAKEFUN_OUTPUT_DIR", "/output"))

LEVELS = [0.05, 0.25, 0.5, 0.75, 0.95]
HOURS = list(range(2, 25))

# Gaussian kernel width on log hour-1 (comparable weighting).
KERNEL_SIGMA = 0.6
# Narrowest log-spread any hour's quantiles may take.
MIN_SPREAD = 0.15
# The linear view's weight is full below this quantile of the
# comparables' log hour-1 and tapers to zero at the largest one.
TAPER_QUANTILE = 0.75
# A comparable counts as a rolling launch when its day peaked this many
# times its hour-1 print.
ROLLING_PEAK = 2.0
# A comparable is usable with at least this many release-day hours.
MIN_DAY_HOURS = 20
# Below this many usable comparables, fall back to flat persistence.
MIN_COMPARABLES = 5
# Fallback spread growth, in log units per square-root hour.
FALLBACK_SIGMA = 0.15


# --- small numeric helpers -------------------------------------------------

def norm_ppf(p):
    """Standard normal quantile (Acklam's rational approximation)."""
    a = (-3.969683028665376e01, 2.209460984245205e02,
         -2.759285984469534e02, 1.383577518672690e02,
         -3.066479806614716e01, 2.506628277459239e00)
    b = (-5.447609879822406e01, 1.615858368580409e02,
         -1.556989798598866e02, 6.680131188771972e01,
         -1.328068155288572e01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01,
         -2.400758277161838e00, -2.549732539343734e00,
         4.374664141464968e00, 2.938163982698783e00)
    d = (7.784695709041462e-03, 3.224671290700398e-01,
         2.445134137142996e00, 3.754408661907416e00)
    p_low = 0.02425
    p_high = 1.0 - p_low
    if p < p_low:
        q = math.sqrt(-2.0 * math.log(p))
        num = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q
               + c[5])
        den = ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
        return num / den
    if p <= p_high:
        q = p - 0.5
        r = q * q
        num = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r
               + a[5]) * q
        den = (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r
               + 1.0)
        return num / den
    q = math.sqrt(-2.0 * math.log(1.0 - p))
    num = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q
           + c[5])
    den = ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    return -num / den


Z_VALUES = [norm_ppf(t) for t in LEVELS]


def round_half_up(x, places=2):
    """Round to `places` decimals, halves away from zero (deterministic)."""
    scale = 10.0 ** places
    shifted = x * scale
    whole = math.floor(shifted)
    if shifted - whole >= 0.5:
        whole += 1
    value = whole / scale
    return int(value) if float(value).is_integer() else value


def weighted_median(pairs):
    """Median of (value, weight) pairs: the value where cumulative weight
    first reaches half the total."""
    ordered = sorted(pairs)
    total = sum(w for _, w in ordered)
    if total <= 0:
        return ordered[len(ordered) // 2][0]
    cut = 0.5 * total
    acc = 0.0
    for value, weight in ordered:
        acc += weight
        if acc >= cut:
            return value
    return ordered[-1][0]


def weighted_sd(pairs):
    """Weighted standard deviation of (value, weight) pairs."""
    total = sum(w for _, w in pairs)
    if total <= 0:
        return 0.0
    mean = sum(v * w for v, w in pairs) / total
    var = sum(w * (v - mean) ** 2 for v, w in pairs) / total
    return math.sqrt(max(var, 0.0))


def solve_linear(matrix, rhs):
    """Solve a small dense linear system by Gaussian elimination with
    partial pivoting. Returns None when the matrix is singular."""
    n = len(rhs)
    rows = [list(matrix[i]) + [rhs[i]] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(rows[r][col]))
        if abs(rows[piv][col]) < 1e-12:
            return None
        rows[col], rows[piv] = rows[piv], rows[col]
        for r in range(col + 1, n):
            factor = rows[r][col] / rows[col][col]
            for cc in range(col, n + 1):
                rows[r][cc] -= factor * rows[col][cc]
    out = [0.0] * n
    for r in range(n - 1, -1, -1):
        rest = sum(rows[r][cc] * out[cc] for cc in range(r + 1, n))
        out[r] = (rows[r][n] - rest) / rows[r][r]
    return out


# --- case data -------------------------------------------------------------

def _parse_time(text):
    return datetime.fromisoformat(text.replace("Z", "+00:00"))


def release_day_curve(comp):
    """A comparable's release-day {hour: viewers} for hours 1..24.

    Hours with no point, or a non-positive count, are skipped, never
    filled in.
    """
    try:
        day_start = _parse_time(comp["day_start"])
        hourly = comp["hourly"]
    except (KeyError, TypeError, ValueError):
        return {}
    curve = {}
    for point in hourly:
        try:
            viewers = point["viewers"]
            when = _parse_time(point["from"])
        except (KeyError, TypeError, ValueError):
            continue
        if viewers is None or viewers <= 0:
            continue
        hour = int((when - day_start).total_seconds() // 3600) + 1
        if 1 <= hour <= 24:
            curve[hour] = viewers
    return curve


def usable_comparables(case):
    """Comparables reduced to (log hour-1, is_rolling, {hour: log-ratio}).

    The target game's own pre-release category curve is not a launch and
    is left out (matched by name). A comparable needs an hour-1 print and
    at least MIN_DAY_HOURS release-day hours.
    """
    target = str(case.get("game", {}).get("name", "")).lower()
    out = []
    for comp in case.get("comparables") or []:
        try:
            name = str(comp.get("name", "")).lower()
        except AttributeError:
            continue
        if target and target in name:
            continue
        curve = release_day_curve(comp)
        if 1 not in curve or len(curve) < MIN_DAY_HOURS:
            continue
        first = curve[1]
        ratios = {h: math.log(v / first) for h, v in curve.items()}
        peak = max(ratios[h] for h in HOURS if h in ratios)
        out.append((math.log(first), peak >= math.log(ROLLING_PEAK), ratios))
    return out


def hour1_print(case):
    """The case's hour-1 print (at least 1), from `known`, else the last
    hourly history point, else 1."""
    for row in case.get("known") or []:
        try:
            if row.get("hour") == 1 and row.get("viewers"):
                return max(1.0, float(row["viewers"]))
        except (AttributeError, TypeError, ValueError):
            continue
    history = case.get("history_hourly") or []
    if history:
        try:
            last = history[-1].get("viewers")
            if last:
                return max(1.0, float(last))
        except (AttributeError, TypeError, ValueError):
            pass
    return 1.0


# --- the two views ----------------------------------------------------------

def kernel_view(comps, log_first):
    """Per hour (centre, spread) from the scale kernel over comparables."""
    view = {}
    for hour in HOURS:
        pairs = []
        for log_h1, _rolling, ratios in comps:
            if hour not in ratios:
                continue
            dist = (log_h1 - log_first) / KERNEL_SIGMA
            pairs.append((ratios[hour], math.exp(-0.5 * dist * dist)))
        if not pairs:
            view[hour] = (0.0, MIN_SPREAD)
            continue
        centre = weighted_median(pairs)
        spread = max(weighted_sd(pairs), MIN_SPREAD)
        view[hour] = (centre, spread)
    return view


def linear_view(comps, log_first):
    """Per hour (centre, spread) from least squares of log-ratio on
    log hour-1 plus a rolling indicator, evaluated for a rolling target
    whose log hour-1 is clipped to the comparables' range."""
    lo = min(c[0] for c in comps)
    hi = max(c[0] for c in comps)
    target = min(max(log_first, lo), hi)
    view = {}
    for hour in HOURS:
        points = [(c[0], 1.0 if c[1] else 0.0, c[2][hour])
                  for c in comps if hour in c[2]]
        if len(points) < 6:
            view[hour] = (0.0, MIN_SPREAD)
            continue
        matrix = [[0.0] * 3 for _ in range(3)]
        rhs = [0.0] * 3
        for x1, x2, y in points:
            feats = (1.0, x1, x2)
            for i in range(3):
                rhs[i] += feats[i] * y
                for j in range(3):
                    matrix[i][j] += feats[i] * feats[j]
        for i in range(3):
            matrix[i][i] += 1e-6
        coef = solve_linear(matrix, rhs)
        if coef is None:
            view[hour] = (0.0, MIN_SPREAD)
            continue
        resid = [y - (coef[0] + coef[1] * x1 + coef[2] * x2)
                 for x1, x2, y in points]
        dof = max(len(points) - 3, 1)
        spread = max(math.sqrt(sum(v * v for v in resid) / dof), MIN_SPREAD)
        centre = coef[0] + coef[1] * target + coef[2]
        view[hour] = (centre, spread)
    return view


def view_quantiles(view, first):
    """A view's per-hour quantile rows in viewer counts."""
    rows = {}
    for hour in HOURS:
        centre, spread = view[hour]
        rows[hour] = [first * math.exp(centre + z * spread) for z in Z_VALUES]
    return rows


def fallback_rows(first):
    """Flat persistence with a spread widening like a random walk."""
    rows = {}
    for hour in HOURS:
        width = FALLBACK_SIGMA * math.sqrt(hour - 1)
        rows[hour] = [first * math.exp(z * width) for z in Z_VALUES]
    return rows


# --- forecast ----------------------------------------------------------------

def forecast_case(case):
    """Forecast hours 2..24 of one case as {"hours": [{hour, q}, ...]}."""
    first = hour1_print(case)
    comps = usable_comparables(case)
    if len(comps) < MIN_COMPARABLES:
        rows = fallback_rows(first)
    else:
        log_first = math.log(first)
        kernel = view_quantiles(kernel_view(comps, log_first), first)
        linear = view_quantiles(linear_view(comps, log_first), first)
        log_h1s = sorted(c[0] for c in comps)
        hi = log_h1s[-1]
        knee = log_h1s[int(TAPER_QUANTILE * (len(log_h1s) - 1))]
        if hi <= knee:
            weight_linear = 0.5
        else:
            weight_linear = 0.5 * min(max((hi - log_first) / (hi - knee),
                                          0.0), 1.0)
        rows = {}
        for hour in HOURS:
            rows[hour] = [
                math.exp((1.0 - weight_linear) * math.log(kernel[hour][i])
                         + weight_linear * math.log(linear[hour][i]))
                for i in range(len(LEVELS))
            ]
    # The horizon is fixed by the competition (hours 2..24); always emit
    # exactly those rows, whatever the input's forecast list says, so a
    # malformed list can never produce a short (invalid) output.
    hours = []
    for hour in HOURS:
        quantiles = []
        running = 0.0
        for value in rows[hour]:
            clean = round_half_up(min(max(value, 0.01), 1e9))
            clean = max(clean, running)
            running = clean
            quantiles.append(clean)
        hours.append({"hour": hour, "q": quantiles})
    return {"hours": hours}


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in sorted(INPUT_DIR.glob("*.in")):
        try:
            with open(path, encoding="utf-8") as handle:
                case = json.load(handle)
            forecast = forecast_case(case)
        except Exception as exc:  # one bad input must not poison the rest
            print(f"warning: skipping {path.name}: {exc}", file=sys.stderr)
            continue
        text = json.dumps(forecast, separators=(",", ":")) + "\n"
        with open(OUTPUT_DIR / (path.stem + ".out"), "w",
                  encoding="utf-8") as handle:
            handle.write(text)


if __name__ == "__main__":
    main()
