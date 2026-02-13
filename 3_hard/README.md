# Hard Test: Regularized Isotonic Regression with Normal Loss

Implement a regularized isotonic regression solver for the Normal/Gaussian loss function.

## Goal

- Implement `NormalLossPiece` class with Constant, Linear, Quadratic coefficients
- Provide `get_smaller_root` and `get_larger_root` implementations
- Use `set_to_min_less_of` for the non-decreasing (isotonic) constraint
- Ideally share a base class with `PoissonLossPieceLog`
- Validate against `isoreg()` (penalty=0) and `Fpop` (non-decreasing changes)
- Write `testthat` unit tests

## Status

TODO
