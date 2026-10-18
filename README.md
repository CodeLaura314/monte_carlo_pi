# Monte Carlo estimation of π

Copyright (C) 2026 lauraoppel

This project contains two small Python programs:

1. **`pi_monte_carlo.py`** — estimates π using a Monte Carlo simulation.
2. **`run_pi_scan.py`** — runs `pi_monte_carlo.py` as a computer experiment in a
   loop, varies the number of samples *n*, analyses the results and generates
   a plot of the error versus *n*.

## The method

Random points `(x, y)` are drawn uniformly from the unit square `[0, 1]²`.
The fraction of points that fall inside the quarter circle of radius 1
(`x² + y² ≤ 1`) approximates the area of the quarter circle, i.e. π/4.
Multiplying that fraction by 4 gives an estimate of π:

```
π̂ = 4 · (points inside the quarter circle) / (total number of samples)
```

The statistical error of the estimator decreases as `1/√n`, which is exactly
what the computer experiment in `run_pi_scan.py` demonstrates.

## `pi_monte_carlo.py`

Runs a single estimation with *n* random samples and prints the estimate
together with its signed error.

```
Usage:
    python3 pi_monte_carlo.py [n]      # default: n = 1000000
```

Example:

```
$ python3 pi_monte_carlo.py 100000
Pi ≈ 3.1344  (nach 100000 Samples, Fehler: -0.007192653589793258)
```

## `run_pi_scan.py`

The computer experiment: calls `pi_monte_carlo.py` in a loop for increasing
sample counts `n` (default: 1 … 10000 in steps of 100), parses the reported
error of each run, prints a table plus a statistical analysis, and plots
`n` vs. error. The plot shows

* the signed error `π̂ − π`,
* the absolute error `|π̂ − π|`,
* a least-squares fit `C/√n`, and
* the theoretical standard deviation `√(4π − π²)/√n` of the estimator,

which together confirm the expected `1/√n` convergence.

Requires **matplotlib** (`pip install matplotlib`).

```
Usage:
    python3 run_pi_scan.py                # scan, print table, save + show plot
    python3 run_pi_scan.py --no-show      # scan, print table, only save plot

Options:
    --start N     first sample count            (default: 1)
    --stop N      last sample count             (default: 10000)
    --step N      step size of the scan         (default: 100)
    --output F    path of the saved plot        (default: pi_error_plot.png)
    --no-show     do not open a plot window
```

The plot is saved as `pi_error_plot.png` next to the script.

## License

GNU General Public License v3.0 or later (see [LICENSE](LICENSE)).

## Requirements

* Python 3
* matplotlib (for `run_pi_scan.py` only)
