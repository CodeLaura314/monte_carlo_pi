#!/usr/bin/env python3
"""Scan n from 1 to 10000 (step 100) by running pi_monte_carlo.py,
parse the reported error and plot n vs. error.

Usage:
    python3 run_pi_scan.py            # scan, print table, save + show plot
    python3 run_pi_scan.py --no-show  # scan, print table, only save plot
"""

# Copyright (C) 2026 lauraoppel
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import argparse
import math
import re
import subprocess
import sys
from pathlib import Path

import matplotlib

# pick the backend BEFORE pyplot is imported (pre-parse --no-show)
if "--no-show" in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROGRAM = Path(__file__).resolve().parent / "pi_monte_carlo.py"
PI_TRUE = 3.141592653589793
ERROR_RE = re.compile(r"Fehler:\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)")


def run_program(n: int) -> float:
    """Run pi_monte_carlo.py with argument n and return the signed error."""
    result = subprocess.run(
        [sys.executable, str(PROGRAM), str(n)],
        capture_output=True,
        text=True,
        check=True,
    )
    match = ERROR_RE.search(result.stdout)
    if not match:
        raise ValueError(f"could not parse output for n={n}: {result.stdout!r}")
    return float(match.group(1))


def scan(start: int, stop: int, step: int) -> tuple[list[int], list[float]]:
    ns, errors = [], []
    for n in range(start, stop + 1, step):
        error = run_program(n)
        ns.append(n)
        errors.append(error)
        print(f"n = {n:6d}   error = {error:+.6e}")
    return ns, errors


def fit_c_over_sqrt_n(ns: list[int], abs_errors: list[float]) -> float:
    """Least-squares fit of |error| = C / sqrt(n)  ->  C = sum(y*x)/sum(x*x)."""
    xs = [1.0 / math.sqrt(n) for n in ns]
    return sum(y * x for x, y in zip(xs, abs_errors)) / sum(x * x for x in xs)


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return cov / (sx * sy)


def analyse(ns: list[int], errors: list[float]) -> None:
    abs_errors = [abs(e) for e in errors]
    mean_abs = sum(abs_errors) / len(abs_errors)
    max_abs = max(abs_errors)
    min_n, max_n = ns[0], ns[-1]
    idx_max = abs_errors.index(max_abs)

    c_fit = fit_c_over_sqrt_n(ns, abs_errors)
    corr = pearson([1.0 / math.sqrt(n) for n in ns], abs_errors)
    # theoretical std of the estimator: sqrt(4*pi - pi^2) / sqrt(n)
    sigma = math.sqrt(4 * PI_TRUE - PI_TRUE**2)

    print("\n=== Analyse ===")
    print(f"Anzahl Messpunkte     : {len(ns)}  (n = {min_n} .. {max_n})")
    print(f"mittlerer |Fehler|    : {mean_abs:.6f}")
    print(f"maximaler |Fehler|    : {max_abs:.6f}  bei n = {ns[idx_max]}")
    print(f"Fit |Fehler| = C/sqrt(n): C = {c_fit:.4f}  (dominiert durch Ausreißer n = {min_n})")
    print(f"theoretische Konstante : sqrt(4*pi - pi^2) = {sigma:.4f}")
    print(f"Korrelation |Fehler| vs. 1/sqrt(n): r = {corr:+.3f}")
    mean_c = sum(a * math.sqrt(n) for a, n in zip(abs_errors, ns)) / len(ns)
    expected_c = sigma * math.sqrt(2 / math.pi)  # E|N(0,sigma)| = sigma*sqrt(2/pi)
    print(f"|Fehler| * sqrt(n)  ->  Mittelwert: {mean_c:.4f}  "
          f"(theoretisch erwartet: {expected_c:.4f})")


def make_plot(ns: list[int], errors: list[float], outfile: Path) -> None:
    abs_errors = [abs(e) for e in errors]
    c_fit = fit_c_over_sqrt_n(ns, abs_errors)
    sigma = math.sqrt(4 * PI_TRUE - PI_TRUE**2)
    curve_n = [n for n in range(max(ns) + 1) if n > 0]
    curve = [c_fit / math.sqrt(n) for n in curve_n]
    sigma_curve = [sigma / math.sqrt(n) for n in curve_n]

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.plot(ns, errors, "o-", markersize=4, alpha=0.65, label="error  $\\hat{\\pi}-\\pi$")
    ax.plot(ns, abs_errors, "s", markersize=4, alpha=0.65, label="absolute error  $|\\hat{\\pi}-\\pi|$")
    ax.plot(curve_n[1:], curve[1:], "r--", linewidth=1.5,
            label=f"fit: {c_fit:.3f}$\\,/\\sqrt{{n}}$")
    ax.plot(curve_n[1:], sigma_curve[1:], "g:", linewidth=1.5,
            label=f"theory: $\\sigma=\\sqrt{{4\\pi-\\pi^2}}/\\sqrt{{n}}$ ({sigma:.3f}$\\,/\\sqrt{{n}}$)")
    ax.plot(curve_n[1:], [-v for v in sigma_curve[1:]], "g:", linewidth=1.5)

    ax.set_xlabel("n (number of samples)")
    ax.set_ylabel("error $\\hat{\\pi} - \\pi$")
    ax.set_title("Monte-Carlo estimation of $\\pi$: error vs. number of samples $n$")
    ax.set_xlim(0, max(ns) * 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    print(f"\nPlot gespeichert: {outfile}")

    # show() on the interactive backend (works after the figure was drawn)
    if matplotlib.get_backend().lower() != "agg":
        plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int, default=1)
    parser.add_argument("--stop", type=int, default=10000)
    parser.add_argument("--step", type=int, default=100)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parent / "pi_error_plot.png")
    parser.add_argument("--no-show", action="store_true", help="do not open a plot window")
    args = parser.parse_args()

    print(f"Scanne n = {args.start} .. {args.stop} mit Schrittweite {args.step}\n")
    ns, errors = scan(args.start, args.stop, args.step)
    analyse(ns, errors)
    make_plot(ns, errors, args.output)


if __name__ == "__main__":
    main()
