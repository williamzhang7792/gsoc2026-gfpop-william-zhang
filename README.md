# GSoC 2026: Time Dependent Constraints in gfpop

**Author:** William (Weidong) Zhang ([GitHub](https://github.com/williamzhang7792))

**Project:** [Time dependent constraints in gfpop](https://github.com/rstats-gsoc/gsoc2026/wiki/Time-dependent-constraints-in-gfpop) -- R Project for Statistical Computing

---

## Test Status

| Test   | Description                                                  | Status    |
|--------|--------------------------------------------------------------|-----------|
| Easy   | `plotModel(graph)` -- visualize model with time-dependent rules | Completed |
| Medium | Optimal partitioning solver with Poisson loss                | TODO      |
| Hard   | Regularized isotonic regression solver with Normal loss      | TODO      |

---

## Easy Test

R function that takes a `gfpop::graph` (with an added `rule` column for the proposed time-dependent constraint API) and draws a state-rule matrix using ggplot2. Rows are states, columns are rules, and arrows show allowed transitions.

See [`1_easy/`](1_easy/) for code and instructions.

![Plot output](1_easy/output_plot.png)

---

## References

- Runge et al., [gfpop: Graph-Constrained Change-point Detection](https://arxiv.org/abs/2002.03646), 2020
- Hocking et al., [LOPART: Labeled Optimal Partitioning](https://doi.org/10.1007/s00180-022-01238-z), 2022
- [gfpop on GitHub](https://github.com/vrunge/gfpop)
