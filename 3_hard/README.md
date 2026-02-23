# Hard Test: Regularized Isotonic Regression (Normal Loss)

A 1-state FPOP solver for regularized isotonic regression with Gaussian (squared-error) loss.

The goal of this implementation is to express isotonic regression within the FPOP framework and make the connection to PAVA explicit.

## Key idea: PAVA ↔ FPOP, same solution in different spaces

PAVA (Pool-Adjacent-Violators Algorithm) enforces monotonicity by pooling adjacent “violators” in **data space**.

The 1-state FPOP formulation enforces the same constraint in **parameter space**.  
The `set_to_min_less_of` operator performs a left-to-right running minimum of the cost function, with transition points determined by quadratic roots. In the end, both approaches produce identical fitted values.

![set_to_min_less_of](set_to_min_less_of.png)

The plot above shows this on a single quadratic cost: the operator traces the input curve (dashed) down to the argmin, then emits a flat line at the minimum cost. This is the geometric mechanism behind PAVA pooling.

## What I implemented

- **`NormalLossPiece`**: quadratic pieces \( f(\mu) = A\mu^2 + B\mu + C \), with closed-form roots.
- **`set_to_min_less_of`**: isotonic operator (running minimum sweep), adapted from gfpop’s “up” operator logic.
- **DP loop (`NormalFPOP`)**: medium-test structure extended with normal loss and the isotonic constraint.
- **Shared base (`LossPiece.h`)**: common interval fields and `clampToInterval`, inherited by both `NormalLossPiece` and the Medium test’s `PoissonLossPieceLog`.

## Design Decisions: OOP vs. Performance

To satisfy the bonus objective, I introduced a `LossPiece` base class to hold shared interval fields and `clampToInterval`.  

However, to preserve empirical performance, `PiecewiseNormalLoss` stores objects using value semantics (`std::list<NormalLossPiece>`) rather than heap-allocated polymorphic pointers. This avoids virtual dispatch overhead inside the inner \( O(N) \) dynamic programming loop while maintaining clean OOP structure.

> If you see a way to simplify the operator logic, I'd love to hear your thoughts. :)

## Tests

`test-NormalFPOP.R` verifies:

- penalty = 0 matches `isoreg()` (includes a 200-trial fuzz test)
- comparisons against `fpop::Fpop()` on monotone-friendly inputs at multiple penalties
- edge cases (n=1, n=2, constant data) and monotonicity invariants

Run:

```bash
Rscript test-NormalFPOP.R