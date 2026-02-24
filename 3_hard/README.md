# Hard Test: Regularized Isotonic Regression (Normal Loss)

This implementation provides a 1-state FPOP solver for regularized isotonic regression with Gaussian (squared-error) loss, expressing isotonic regression within the FPOP framework and making the connection to PAVA explicit.

## Key idea: PAVA and FPOP solve the same problem in different spaces

PAVA (Pool-Adjacent-Violators) enforces monotonicity by pooling adjacent violators in **data space** — it merges neighboring segments whenever the mean decreases.

The 1-state FPOP formulation enforces the same constraint in **parameter space**. The `set_to_min_less_of` operator sweeps a running minimum over the piecewise quadratic cost from left to right: it traces the curve downhill, then emits a flat line once the minimum is reached. Transition points between "follow the curve" and "flatline" are determined by quadratic roots.

With penalty = 0, both approaches produce identical fitted values. Increasing the penalty regularizes the number of segments, which PAVA alone cannot do.

![set_to_min_less_of](set_to_min_less_of.png)

## Approach

I started from my Medium test solver and made two changes:

1. **Swapped the loss**: replaced `PoissonLossPieceLog` (exp/log pieces) with `NormalLossPiece` (quadratic pieces, `f(μ) = Aμ² + Bμ + C`). Quadratic pieces have closed-form roots, so no Newton iteration is needed.
2. **Swapped the constraint operator**: replaced `set_to_unconstrained_min_of` (global min → flat piece) with `set_to_min_less_of` (running min sweep, adapted from gfpop's "up" operator). This is the isotonic constraint — it only allows the optimal mean to stay the same or increase.

Everything else carries over: the FPOP pruning rule, `set_to_min_env_of` for combining cost functions, `Minimize` for extracting the optimal mean, and the backtracking loop.

## Design decision: OOP vs. performance

To satisfy the bonus objective, both `NormalLossPiece` and `PoissonLossPieceLog` inherit from a shared `LossPiece` base class (interval fields + `clampToInterval`).

However, `PiecewiseNormalLoss` stores pieces by value (`std::list<NormalLossPiece>`) rather than through polymorphic pointers. This avoids virtual dispatch overhead inside the O(N) inner loop while still keeping a clean class hierarchy.

> If you see a way to simplify the operator logic, I'd love to hear your thoughts. :)

## Quickstart

Dependencies:

```r
install.packages(c("Rcpp", "testthat"))
# optional, for comparison tests:
install.packages("fpop")
```

Run:

```bash
cd 3_hard
Rscript test-NormalFPOP.R
```

Compiles `NormalFPOP.cpp` via `Rcpp::sourceCpp` and runs the full test suite.

## Validation

- **penalty = 0**: fitted values match `isoreg()` exactly, verified on hand-picked cases and a 200-trial fuzz test (N = 10, mixed distributions).
- **penalty > 0**: on non-decreasing inputs (where the isotonic constraint is not active), results match `fpop::Fpop()` across penalties 1, 5, 10, 50.
- **Invariants**: segment means are always non-decreasing; high penalty collapses to 1 segment at the overall mean.
- **Edge cases**: n = 1, n = 2, constant data, already-sorted data, strictly decreasing data.

## Files

| File | Description |
|------|-------------|
| `NormalFPOP.cpp` | Rcpp solver: `NormalLossPiece`, `PiecewiseNormalLoss`, `set_to_min_less_of`, DP loop |
| `LossPiece.h` | Shared base class, also used by Medium test's `PoissonLossPieceLog` |
| `test-NormalFPOP.R` | `testthat` suite: isoreg comparison, fpop comparison, fuzz test, edge cases |
| `figure_set_to_min_less_of.R` | Generates the operator visualization above |
| `set_to_min_less_of.png` | Generated figure |
| `reference/` | Study copies of gfpop source I read while working (gitignored) |

## References

- Runge et al. (2020), *gfpop: Graph-Constrained Change-point Detection*  
  https://arxiv.org/abs/2002.03646
- Maidstone et al. (2016), *On Optimal Multiple Changepoint Algorithms for Large Data*  
  https://doi.org/10.1007/s11222-016-9636-3
- PeakSegOptimal (GitHub)  
  https://github.com/tdhock/PeakSegOptimal
- Barlow et al. (1972), *Statistical Inference under Order Restrictions* — isotonic regression reference
- `isoreg` (R base)  
  https://stat.ethz.ch/R-manual/R-devel/library/stats/html/isoreg.html
