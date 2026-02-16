# Easy Test: `plotModel(graph)`

Visualize a changepoint model with time-dependent constraints as a state-rule matrix.

## How to Run

Install dependencies (gfpop is only on GitHub, not CRAN):

```r
install.packages("remotes")
remotes::install_github("vrunge/gfpop")
install.packages("ggplot2")
```

Then run:

```bash
Rscript plotModel.R
```

This generates `output_plot.png` in the current directory.

## What It Does

- `create_LOPART_graph()` builds a graph using `gfpop::Edge()` calls matching the wiki example, then adds a `rule` column manually (since gfpop doesn't support it yet).
- `plotModel(graph)` takes that graph and draws a ggplot2 visualization: rows = states, columns = rules, arrows = allowed transitions. Self-loops are drawn as curved arcs above the node.

## Output

![output_plot.png](output_plot.png)
