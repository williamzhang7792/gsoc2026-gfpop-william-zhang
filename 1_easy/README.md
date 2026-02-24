# Easy Test: `plotModel(graph)`

This script visualizes a changepoint model with time-dependent constraints as a state–rule transition matrix.

## Approach

The proposed gfpop extension introduces a `rule` column to the graph, making transitions time-dependent. `plotModel` takes a graph with this column and renders it as a matrix: rows = states, columns = rules, arrows = allowed transitions per time interval.

I built the LOPART graph from the GSoC wiki example using real `gfpop::Edge()` calls, then added the `rule` column manually (since gfpop doesn't support it yet). This shows what the extended API would look like in practice.

## Quickstart

Dependencies (gfpop is only on GitHub, not CRAN):

```r
install.packages("remotes")
remotes::install_github("vrunge/gfpop")
install.packages("ggplot2")
```

Run:

```bash
cd 1_easy
Rscript plotModel.R
```

Produces `output_plot.png`.

## What it does

- `create_LOPART_graph()` builds an 8-edge graph via `gfpop::Edge()` and attaches a `rule` column with semantic labels (unlabeled, in positive label, etc.)
- `plotModel(graph)` maps states to rows and rules to columns, then draws arrows for each allowed transition. Self-loops render as curved arcs above the node; penalty labels appear on edges where applicable.
- Nodes are drawn last (white-filled labels) so arrowheads terminate cleanly underneath.

## Workaround

gfpop doesn't have a `rule` column yet — that's the proposed extension. I add it manually after graph construction and store rule labels as an attribute. This keeps the demo functional without modifying the package.

## Files

- `plotModel.R` — entrypoint; defines `create_LOPART_graph()`, `plotModel()`, and runs both
- `requirements.txt` — R package dependencies
- `output_plot.png` — generated output

## Output

![output_plot.png](output_plot.png)

## References

- Runge et al. (2020), *gfpop: Graph-Constrained Change-point Detection*  
  https://arxiv.org/abs/2002.03646
- Hocking et al. (2022), *LOPART: Labeled Optimal Partitioning*  
  https://doi.org/10.1007/s00180-022-01238-z
- gfpop (GitHub)  
  https://github.com/vrunge/gfpop
