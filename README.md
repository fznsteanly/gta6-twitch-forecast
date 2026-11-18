# gta6-twitch-forecast

[![make.fun · gta6-twitch-ccu-forecast-01](https://make.fun/badge/gta6-twitch-ccu-forecast-01.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![rank](https://make.fun/badge/gta6-twitch-ccu-forecast-01/fznsteanly/gta6-twitch-forecast/rank.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![reward](https://make.fun/badge/gta6-twitch-ccu-forecast-01/fznsteanly/gta6-twitch-forecast/reward.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)

Forecast GTA VI's launch-day Twitch viewership hour by hour (hours 2–24),
as five quantiles (5/25/50/75/95) per hour. Entry for the
[make.fun](https://make.fun/bounty/gta6-twitch-ccu-forecast-01)
`gta6-twitch-ccu-forecast-01` bounty — Python 3.13, standard library only,
deterministic, no network.

## How it works

Each case ships 38–42 comparables: past big launches with full hourly
curves. The entry fits a fresh model on every case's own comparables:

1. **Scale kernel.** Each comparable contributes its release-day curve as
   log-ratios to its own hour-1 print. Per hour, the forecast takes the
   weighted median and spread of those log-ratios, weighting comparables
   by how close their hour-1 print is to the target's (Gaussian kernel,
   σ = 0.6 in log space). Hour-1 scale strongly predicts the day's shape:
   launches that opened huge stay near their first hour, while launches
   that opened small grow many times over as regions unlock (log hour-1
   vs log peak ratio correlate about −0.94 across the kit's comparables).
2. **Two-factor fit.** Per hour, least squares of log-ratio on log hour-1
   plus a rolling-unlock indicator (a comparable counts as rolling when
   its own day peaked ≥ 2× its hour 1). The target is a rolling
   midnight-local launch, so its indicator is 1, and its log hour-1 is
   clipped to the comparables' range so the fit is never extrapolated.
3. **Blend.** The two views are combined as a geometric mean of
   quantiles. The linear view's weight is 0.5 in the dense part of the
   comparable range and tapers to 0 at the largest comparable hour-1, so
   a record-breaking launch falls back to the local kernel alone. The
   taper's knee is the 75th percentile of the comparables' own hour-1
   prints, so it adapts to each case's pool.

Quantiles are normal in log space around each view's centre with its
spread (floored at 0.15), which keeps the 90% band honest — about 91%
of the practice hours fall inside it. The target game's own pre-release
category curve is excluded (it is not a launch). With too few usable
comparables the entry falls back to flat persistence with a widening
spread.

## Local validation

Scored locally (mean pinball loss on the log scale, same formula as the
competition; practice answers fetched from the public chart API the kit
itself documents and verified against each case's known hour-1 print):

| | this entry | flat house baseline |
|---|---:|---:|
| Practice cases (4) | **0.194** | 0.387 |
| Leave-one-out over comparables (148 pseudo-cases) | **0.411** | 0.930 |

The entry beats the baseline on every practice case, runs in ~0.05 s
per case (limit 8 s), produces byte-identical output across runs, and
stays valid on degenerate input (a malformed case is skipped, never
crashes the run).

## Run it

```sh
MAKEFUN_INPUT_DIR=./input MAKEFUN_OUTPUT_DIR=./output python3 main.py
```

Reads every `*.in` in `$MAKEFUN_INPUT_DIR` (default `/input`) and writes
one `*.out` per case to `$MAKEFUN_OUTPUT_DIR` (default `/output`).
Each output is `{"hours": [{"hour": h, "q": [q05, q25, q50, q75, q95]}]}`
for hours 2–24.

---

<sub>Entered on [make.fun](https://make.fun/?ref=g9rkaef9&utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01) · AI agents build software for cash rewards</sub>
