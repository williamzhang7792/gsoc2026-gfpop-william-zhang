# GSoC 2026 Qualification Tests: Time Dependent Constraints in gfpop

**Author:** William (Weidong) Zhang -- [GitHub](https://github.com/williamzhang7792)

**Project:** [Time dependent constraints in gfpop](https://github.com/rstats-gsoc/gsoc2026/wiki/Time-dependent-constraints-in-gfpop) -- R Project for Statistical Computing

**Repository:** [gsoc2026-gfpop-william-zhang](https://github.com/williamzhang7792/gsoc2026-gfpop-william-zhang)

---

## Test Status

| Test   | Description                                                                 | Status         |
|--------|-----------------------------------------------------------------------------|----------------|
| Easy   | `plotModel(graph)` -- visualize changepoint model with time-dependent rules | Completed |
| Medium | Optimal partitioning solver with Poisson loss (unconstrained)               | TODO           |
| Hard   | Regularized isotonic regression solver with Normal loss                     | TODO           |

---

## Easy Test: `plotModel`

An R function that takes a `gfpop::graph`-style dataframe (with an added `rule` column for time-dependent constraints) and visualizes it as a state-rule matrix using `ggplot2`.

- **Rows (Y-axis):** States (e.g., `"normal"`, `"noChange"`)
- **Columns (X-axis):** Rules (integer IDs specifying which edges are active at each time point)
- **Edges:** Directed arrows between nodes; self-loops rendered as curved arcs
- **Labels:** Each edge is annotated with its cost type and penalty

### Output

![Changepoint Model Plot](1_easy/output_plot.png)

### Files

| File                                       | Description                         |
|--------------------------------------------|-------------------------------------|
| [`1_easy/plotModel.R`](1_easy/plotModel.R)     | Main script with `plotModel()` function and mock data |
| [`1_easy/output_plot.png`](1_easy/output_plot.png) | Generated plot                  |
| [`1_easy/requirements.txt`](1_easy/requirements.txt) | R package dependencies        |

See [`1_easy/README.md`](1_easy/README.md) for instructions on running the script.

---

## Medium Test

*TODO*

---

## Hard Test

*TODO*

---

## Project Background

The [gfpop](https://github.com/vrunge/gfpop) package implements Graph-Constrained Functional Pruning Optimal Partitioning for changepoint detection. Currently it only supports constraints that are valid for **all** time points. The goal of this GSoC project is to extend gfpop to allow **time-dependent constraints**, enabling fast optimal solvers for models such as Labeled Optimal Partitioning (LOPART).

The key idea is to add a `rule` argument so that different subsets of graph edges are used at different data points, controlled by an integer vector of rule IDs.

## References

- Runge et al., [gfpop: An R Package for Univariate Graph-Constrained Change-point Detection](https://arxiv.org/abs/2002.03646), 2020
- Hocking et al., [LOPART: Labeled Optimal Partitioning](https://doi.org/10.1007/s00180-022-01238-z), 2022
- [gfpop on CRAN](https://cran.r-project.org/package=gfpop) | [gfpop on GitHub](https://github.com/vrunge/gfpop)
