# Hard Test: Regularized Isotonic Regression (Normal Loss)

A 1-state FPOP solver for regularized isotonic regression with Gaussian loss,
adapted from the Medium test (`PoissonFPOP.cpp`).

## What's implemented

- **`NormalLossPiece`** class with `Quadratic`, `Linear`, `Constant` coefficients
  for `f(mu) = Q*mu^2 + L*mu + C`. Root-finding via the quadratic formula (`sqrt`
  from `<math.h>`).
- **`set_to_min_less_of`** operator enforcing the non-decreasing (isotonic)
  constraint, adapted from gfpop's `operatorUp` / `intervalMinLessUp` /
  `pastePieceUp`.
- **Shared base class** (`LossPiece` in `LossPiece.h`): both `NormalLossPiece`
  and `PoissonLossPieceLog` (in `PoissonFPOP.cpp`) inherit interval fields and
  `clampToInterval`.
- **`NormalFPOP(data, penalty)`** Rcpp export returning segment means, ends,
  breakpoints, cost, and interval counts.

## Tests

`test-NormalFPOP.R` covers:
- penalty=0 vs `isoreg()` on sorted, decreasing, constant, and random data
- non-decreasing data vs `fpop::Fpop()` at several penalties
- edge cases (n=1, n=2, high penalty, isotonic invariant)
- 200-trial random fuzz (normal, uniform-integer, large-magnitude distributions)

Run:
```
Rscript test-NormalFPOP.R
```

## Status

Complete.
