# Medium Test: Unconstrained FPOP for Poisson Loss

This implementation applies the FPOP algorithm (Maidstone et al. 2016) to optimal partitioning with Poisson loss and no constraints, adapted from PeakSegOptimal's constrained solver.

For a step-by-step walkthrough of the algorithm with a worked example, see **[algo.pdf](algo.pdf)** (compiled from `algo.tex`).

## Approach

PeakSegOptimal solves a 2-state (up/down) constrained problem using `set_to_min_less_of` / `set_to_min_more_of` to enforce monotonicity between adjacent segments. I stripped out the constraint layer to get a 1-state unconstrained solver:

1. Replaced both monotonicity operators with `set_to_unconstrained_min_of` — this just finds the global minimum and emits a single flat piece.
2. Collapsed the cost array from N × 2 (up/down states) to N × 1.
3. Kept everything else intact: the piecewise function infrastructure (`set_to_min_env_of`, `Minimize`, Newton root-finding), the FPOP pruning rule, and the backtracking logic.

The key insight is that unconstrained FPOP is a strict simplification of the constrained case — removing constraints means fewer pieces to track per step, so the pruning is even more aggressive.

## Quickstart

Dependencies:

```r
install.packages(c("Rcpp", "testthat"))
remotes::install_github("cran/Segmentor3IsBack")
```

Run:

```bash
cd 2_medium
Rscript PoissonFPOP.R        # comparison with Segmentor3IsBack
Rscript test-PoissonFPOP.R   # testthat suite
```

`PoissonFPOP.R` compiles the solver via `Rcpp::sourceCpp` and prints a penalty-by-penalty comparison. `test-PoissonFPOP.R` runs 9 test cases.

To recompile the algorithm walkthrough:

```bash
tectonic algo.tex    # or: pdflatex algo.tex
```

## Implementation note

The DP loop divides by cumulative weight at each step (`cost->multiply(1.0/cum_weight_i)`). This normalization is inherited from PeakSegOptimal's convention for the Poisson loss parameterization in log-mean space; the Hard test's Normal loss solver does not need it because the quadratic cost is expressed directly in mean space.

## Validation

Breakpoints and segment means match `Segmentor3IsBack::Segmentor(model=1)` exactly across all test cases: 3-segment Poisson data, single changepoint, constant rate, 5-segment, data with zeros, penalty=0, n=2, and bad-input handling.

## Files

| File | Description |
|------|-------------|
| `PoissonFPOP.cpp` | Self-contained Rcpp solver (~760 lines); `PoissonLossPieceLog` inherits from `LossPiece` |
| `LossPiece.h` | Shared base class (interval fields + `clampToInterval`), also used by Hard test |
| `PoissonFPOP.R` | Compiles solver, runs side-by-side comparison with Segmentor3IsBack |
| `test-PoissonFPOP.R` | `testthat` suite: 9 test cases |
| `algo.tex` / `algo.pdf` | Algorithm walkthrough with worked example (Beamer slides) |
| `reference/` | Study copies of PeakSegOptimal source I read while working |

## Output

![FPOP Pruning Trace](pruning_demo.png)

## References

- Maidstone et al. (2016), *On Optimal Multiple Changepoint Algorithms for Large Data*  
  https://doi.org/10.1007/s11222-016-9636-3
- PeakSegOptimal (GitHub)  
  https://github.com/tdhock/PeakSegOptimal
- Segmentor3IsBack (CRAN)  
  https://cran.r-project.org/package=Segmentor3IsBack
