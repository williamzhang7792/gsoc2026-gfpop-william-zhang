# Easy Test: `plotModel(graph)`

Visualize a changepoint model with time-dependent constraints as a state-rule matrix.

## Requirements

- **R** (>= 4.0)
- **gfpop** (>= 1.1.1) -- install from GitHub: `remotes::install_github("vrunge/gfpop")`
- **ggplot2** (>= 3.0)

Install dependencies from R:

```r
install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("vrunge/gfpop")
install.packages("ggplot2", repos = "https://cloud.r-project.org")
```

## How to Run

From this directory:

```bash
Rscript plotModel.R
```

This will generate `output_plot.png` in the current directory.

Alternatively, from R or RStudio:

```r
source("plotModel.R")
run_main()
```

## What It Does

1. **`create_LOPART_graph()`** -- builds a real `gfpop::graph` using `gfpop::Edge()` calls, then manually adds a `rule` column (not yet supported by gfpop) to demonstrate the proposed time-dependent constraint API (4 rules, 2 states).
2. **`plotModel(graph)`** -- takes the graph dataframe and produces a `ggplot2` visualization:
   - Nodes are arranged in a matrix (rows = states, columns = rules).
   - Directed arrows show allowed transitions.
   - Self-loops are drawn as curved arcs with alternating direction.
   - Edge labels show the cost type and penalty.

## Output

![output_plot.png](output_plot.png)

## Reproducibility

Generated with:

- R 4.5.1 (2025-06-13)
- ggplot2 4.0.2
- Platform: aarch64-apple-darwin20 (macOS ARM64)
