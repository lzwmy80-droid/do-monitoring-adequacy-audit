# Nonnegative state constraint (`x(t) >= 0`)

Status: `implemented_in_current_runtime_and_checked_against_small_n_oracle`.

## Decision

For dissolved oxygen, add the physical state space

```text
F_+(S) = {x: x is globally L-Lipschitz, x(t_i)=y_i for i in S, x(t)>=0}.
```

Here `S` is one feasible retained-index set after at most `k` replacements.
This is not a cosmetic truncation of a reported interval. It changes the
feasible trajectory class and therefore must be stated in the theorem,
algorithm, and empirical estimand.

The change is technically small. For nonnegative retained reports, the exact
fixed-`S` constrained envelopes are

```text
E_S^+(t) = max(E_S(t), 0)
U_S^+(t) = U_S(t),
```

where `E_S` and `U_S` are the existing McShane--Whitney lower and upper
envelopes.

## Why direct clipping is exact and sharp

Every constrained feasible path already obeys `x >= E_S` and `x >= 0`, hence
`x >= max(E_S,0)`. The scalar projection `z -> max(z,0)` is 1-Lipschitz, so
`max(E_S,0)` remains globally `L`-Lipschitz. It agrees with every retained
anchor because `y_i >= 0`. Also `U_S >= 0`, because it is the minimum of cones
`y_i + L|t-t_i|`, and `max(E_S,0) <= U_S`. Thus the clipped lower envelope and
the unchanged upper envelope are both feasible and attain the pointwise
extrema.

For a positive environmental threshold `H>0`, clipping leaves the strict
sublevel set unchanged:

```text
{max(E_S,0) < H} = {E_S < H}.
```

For cumulative deficit,

```text
(H-max(E_S,0))_+ = (H-E_S)_+ - (0-E_S)_+.
```

Therefore, for a fixed feasible `S`:

| Endpoint | Effect of `x>=0` when `H>0` |
|---|---|
| occupation lower, attained by `U_S` | unchanged |
| occupation upper, attained by `E_S` | unchanged |
| deficit lower, attained by `U_S` | unchanged |
| deficit upper, attained by `E_S` | replace `E_S` by `max(E_S,0)` |

If every report is nonnegative, the admissible retained subsets are unchanged.
After optimizing over the replacement union, the first three robust endpoints
therefore remain unchanged. The deficit-upper retained/deleted witness may
change because its edge costs change, but it is still attained by the clipped
lower envelope for its selected subset.

If negative reports are allowed as inputs, they cannot be retained. They are
mandatory replacements, so `k_min(L)` and potentially all four optimized
endpoints can change. The frozen MARACOOS preprocessing already keeps DO in
`[0,50] mg/L`, so this issue does not affect that dataset.

## Immediate scientific consequence

For `H>0`,

```text
0 <= D_H(x) <= H*T
0 <= D_upper - D_lower <= H*T.
```

Hence a reported `deficit width/(H*T)` above one is impossible under the
nonnegative physical model. Values above one in the current empirical tables
are mathematically valid for the unconstrained real-valued signal model, but
they are driven partly by negative between-sample trajectories and should not
be presented as the final dissolved-oxygen result.

When `H<=0`, both occupation below `H` and threshold deficit are identically
zero under `x>=0`.

## No-anchor branch

If `k>=n` and every report may be replaced, the unconstrained runtime currently
returns an unattained infinite deficit upper endpoint. Under `x>=0` this
becomes

```text
occupation = [0,T] and deficit = [0,H*T]  when H>0,
occupation = [0,0] and deficit = [0,0]    when H<=0.
```

The upper endpoints are attained by the constant zero path. Thus every
nonnegative four-endpoint witness is finite and attained.

## Current production implementation

`robust_exact_current.py` now exposes the optional keyword
`state_lower_bound`; the unconstrained behavior remains the default and the DO
analysis passes `state_lower_bound=0.0`.

The current implementation:

1. excludes reports below the floor from every retained-subsequence state and
   recomputes the floor-aware `k_min`;
2. clips lower-envelope affine pieces at the floor before occupation/deficit
   integration;
3. reoptimizes every endpoint with the floor-aware additive edge costs, rather
   than clipping only the old unconstrained witness;
4. preserves `O(k n^2)` time and `O(k n)` dynamic-programming memory;
5. returns finite attained endpoints when all reports are replaceable; and
6. labels the result mode with the numeric floor, for example
   `exact_k_replacement_contamination_with_state_floor=0`.

Interval-observation support is outside this exact-inlier runtime and is not
silently implied by the new option.

## Required regression gates

The independent exhaustive oracle is
`nonnegative_state_constraint_prototype.py`; it is not production code.
`test_nonnegative_state_constraint_prototype.py` covers:

1. a two-anchor V-shaped lower envelope where only deficit upper changes
   (`2.0 -> 1.5`);
2. exact equality with the old runtime when the lower envelope never crosses
   zero;
3. a finite, attained `H*T` upper deficit when all reports are replaceable;
4. mandatory deletion of a negative exact report;
5. collapse of both functionals when `H<=0`;
6. equality between explicit clipping and the hinge-difference identity;
7. random small replacement cases in which all four production endpoints
   match the independent exhaustive floor-constrained oracle and constrained
   deficit width never exceeds `H*T`.

The state-floor result is a new estimand. It must be written to separate
outputs rather than overwriting the existing unconstrained audit tables.
