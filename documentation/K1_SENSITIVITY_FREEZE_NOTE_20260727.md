# Fixed-\(k=1\) sensitivity analysis freeze note

Status: analysis specification frozen before the first \(k=1\) computation.  
This is an internal analysis freeze, not a prospective registration.

## Purpose

The seven-day Table 2 experiment isolates sampling uncertainty at \(k=0\).
This sensitivity analysis asks how much the sharp endpoint hull expands under
the smallest positive, externally fixed replacement budget, \(k=1\). The
budget is not estimated from a case and is not set equal to
\(k_{\min}(L)\).

## Frozen scope

- Use exactly the 240 seven-day station-window-cadence-phase cases already
  classified as feasible and reference-\(L\)-compatible in the frozen
  \(k=0\) analysis.
- Preserve the holdout split, 15-minute operational reference, all cadence
  phases, \(H=5\) mg/L, calibration-only q99.9 value of \(L\), and physical
  state floor \(B=0\) mg/L.
- Change only the maximum replacement budget from \(k=0\) to \(k=1\).
- Do not retune \(L\), select stations, or exclude cases after observing the
  \(k=1\) results.

## Primary summaries

For each cadence (0.5, 1, 2, and 4 h), report:

1. the number of analyzed cases;
2. median normalized occupation width at \(k=0\) and \(k=1\);
3. median normalized nonnegative deficit width at \(k=0\) and \(k=1\);
4. median within-case width increase for each functional;
5. the number of cases in which at least one endpoint certificate uses the
   permitted deletion;
6. reference-containment and nesting-check failures.

Normalization remains \(T\) for occupation and \(HT\) for cumulative deficit.
The operational reference is used only for deterministic containment checks,
not as statistical truth or for tuning.

## Interpretation rule

The \(k=1\) result is a fixed-budget sensitivity scenario. It is not an
estimated corruption prevalence and does not establish that any deleted report
is erroneous.
