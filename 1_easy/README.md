# Easy Test: `plotModel(graph)`

Visualize a changepoint model with time-dependent constraints as a state-rule matrix.

## Requirements

- **R** (>= 4.0)
- **ggplot2** (>= 3.0)

Install the dependency from R:

```r
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

1. **`create_mock_graph()`** -- builds a dataframe mimicking `gfpop::graph()` output with an added `rule` column, encoding the LOPART model (4 rules, 2 states).
2. **`plotModel(graph)`** -- takes the graph dataframe and produces a `ggplot2` visualization:
   - Nodes are arranged in a matrix (rows = states, columns = rules).
   - Directed arrows show allowed transitions.
   - Self-loops are drawn as curved arcs with alternating direction.
   - Edge labels show the cost type and penalty.

## Output

![output_plot.png](output_plot.png)
