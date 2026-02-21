# GSoC 2026 — Time-Dependent Constraints in gfpop

This repository contains my attempts at the Easy, Medium, and Hard tests for the gfpop project.

---

## Progress

- ✅ Easy — `plotModel(graph)` visualization
- ✅ Medium — Unconstrained FPOP (Poisson loss)
- 🚧 Hard — Regularized isotonic regression (Normal loss, in progress)

---

## Easy Test

Implemented an R function that takes a `gfpop::graph` (extended with a `rule` column) and visualizes the state–rule transition matrix using ggplot2.

Rows represent states, columns represent rules, and arrows indicate allowed transitions.

See [`1_easy/`](1_easy/) for code and instructions.

![Plot output](1_easy/output_plot.png)

---

## Medium Test

Implemented an unconstrained FPOP solver for optimal partitioning with Poisson loss. I adapted the `PeakSegOptimal` architecture from a 2-state constrained model to a 1-state unconstrained dynamic program. The implementation was validated against `Segmentor3IsBack`, matching its output while achieving empirical $O(N)$ performance.

See [`2_medium/`](2_medium/) for code and usage instructions. I also prepared a short slide deck ([`algo.pdf`](2_medium/algo.pdf), compiled from `algo.tex`) outlining my current understanding of the FPOP pruning mechanics.

Feedback or corrections are very welcome — I’m always looking to get a better understanding. :)

![FPOP Pruning Trace](2_medium/pruning_demo.png)

---

## References

- Runge et al., [gfpop: Graph-Constrained Change-point Detection](https://arxiv.org/abs/2002.03646), 2020
- Maidstone et al., [On Optimal Multiple Changepoint Algorithms for Large Data](https://doi.org/10.1007/s11222-016-9636-3), 2017
- Hocking et al., [LOPART: Labeled Optimal Partitioning](https://doi.org/10.1007/s00180-022-01238-z), 2022
- [PeakSegOptimal on GitHub](https://github.com/tdhock/PeakSegOptimal)
- [Segmentor3IsBack on CRAN](https://cran.r-project.org/package=Segmentor3IsBack)
- [gfpop on GitHub](https://github.com/vrunge/gfpop)
