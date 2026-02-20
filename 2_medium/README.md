# Medium Test: Unconstrained FPOP for Poisson Loss

Implements the FPOP algorithm (Maidstone et al. 2016) for optimal partitioning with Poisson loss and no constraints, by adapting the constrained (up-down) solver from [PeakSegOptimal](https://github.com/tdhock/PeakSegOptimal).

For a step-by-step walkthrough of the algorithm with a worked example, see **[algo.pdf](algo.pdf)** (compiled from `algo.tex`).

## Approach

The constrained solver uses two states (up/down) with `set_to_min_less_of` / `set_to_min_more_of` to enforce monotonicity between adjacent segment means. The unconstrained solver:

1. Replaces both with `set_to_unconstrained_min_of` (global minimum → flat piece)
2. Simplifies the cost array from N x 2 to N x 1
3. Reuses all piecewise function infrastructure (`set_to_min_env_of`, `Minimize`, root-finding, etc.)

## Files

| File | Description |
|------|-------------|
| `PoissonFPOP.cpp` | Self-contained Rcpp solver (~760 lines, mostly reused from `funPieceListLog`) |
| `PoissonFPOP.R` | Compiles the solver and compares output with Segmentor3IsBack |
| `test-PoissonFPOP.R` | `testthat` suite: 8 test cases, 54 assertions |
| `algo.tex` / `algo.pdf` | Algorithm walkthrough with worked example (Beamer) |
| `reference/` | Study copies of PeakSegOptimal source |

## How to run

```r
install.packages(c("Rcpp", "testthat"))
remotes::install_github("cran/Segmentor3IsBack")

setwd("2_medium")
source("PoissonFPOP.R")      # comparison script
source("test-PoissonFPOP.R") # tests
```

```bash
cd 2_medium
Rscript PoissonFPOP.R
Rscript test-PoissonFPOP.R
```

To recompile the algorithm walkthrough:

```bash
tectonic algo.tex    # or: pdflatex algo.tex
```

## Results

All tests pass. Breakpoints and segment means match `Segmentor3IsBack::Segmentor(model=1)` exactly across 8 test cases (3-segment data, single changepoint, constant rate, 5-segment, zeros, penalty=0, n=2, bad input).

## Reference

- Original solver: [PeakSegFPOPLog.cpp](https://github.com/tdhock/PeakSegOptimal/blob/master/src/PeakSegFPOPLog.cpp)
- Algorithm: Maidstone et al. 2016, "On Optimal Multiple Changepoint Algorithms for Large Data"
- Ground truth: [Segmentor3IsBack](https://cran.r-project.org/package=Segmentor3IsBack)
