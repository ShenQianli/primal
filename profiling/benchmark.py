"""Performance benchmark for pyprimal solvers.

Measures CPU time of Dantzig selector and compressed sensing across
increasing dimensions, mirroring the R benchmark (benchmark.R).

Usage:
    python benchmark.py              # default dimensions
    python benchmark.py --large      # larger dimensions (slower, matches R benchmark range)
"""

import argparse
import sys
import time

import numpy as np

from pyprimal import dantzig_solver, compressed_sensing_solver


def generate_data(n, d, sparsity=0.03, seed=None):
    """Generate sparse regression data (matches R dantzig.generator)."""
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, d))
    s = max(1, int(d * sparsity))
    beta = np.zeros(d)
    beta[:s] = rng.uniform(-1, 1, s)
    y = X @ beta + rng.standard_normal(n)
    return X, y


parser = argparse.ArgumentParser(description="Benchmark pyprimal solvers")
parser.add_argument(
    "--large", action="store_true",
    help="Use large dimension range (d up to 7000, much slower)",
)
args = parser.parse_args()

if args.large:
    dims_dantzig = list(range(1000, 8000, 1000))  # matches R benchmark
    dims_cs = list(range(200, 1100, 100))
    trial_n = 5
else:
    dims_dantzig = list(range(100, 600, 100))
    dims_cs = list(range(100, 600, 100))
    trial_n = 3

n = 200


# ── Dantzig selector benchmark ──────────────────────────────────────

print("************** Benchmark: Dantzig selector ****************")

times_dantzig = []

for d in dims_dantzig:
    lam = 2 * np.sqrt(np.log(d) / n)
    X, y = generate_data(n, d, sparsity=0.03, seed=42)

    t0 = time.perf_counter()
    dantzig_solver(X, y, max_it=10000, lambda_threshold=lam)
    elapsed = time.perf_counter() - t0
    times_dantzig.append(elapsed)

    print(f"  d={d:5d}  time={elapsed:.4f}s")
    sys.stdout.flush()


# ── Compressed sensing benchmark ─────────────────────────────────────

print("\n************** Benchmark: Compressed sensing ***************")

times_cs = []

for d in dims_cs:
    lam = 2 * np.sqrt(np.log(d) / n)
    total = 0.0

    for trial in range(trial_n):
        X, y = generate_data(n, d, sparsity=0.1, seed=42 + trial)

        t0 = time.perf_counter()
        compressed_sensing_solver(X, y, max_it=10000, lambda_threshold=lam)
        total += time.perf_counter() - t0

    avg = total / trial_n
    times_cs.append(avg)
    print(f"  d={d:4d}  avg_time={avg:.4f}s  ({trial_n} trials)")
    sys.stdout.flush()


# ── Plot results ─────────────────────────────────────────────────────

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].plot(dims_dantzig, times_dantzig, "b-o", linewidth=2, label="pyprimal")
    axes[0].set_title("Dantzig Selector", fontsize=14)
    axes[0].set_xlabel("Dimension", fontsize=12)
    axes[0].set_ylabel("CPU Time (s)", fontsize=12)
    axes[0].legend(fontsize=12)

    axes[1].plot(dims_cs, times_cs, "b-o", linewidth=2, label="pyprimal")
    axes[1].set_title("Compressed Sensing", fontsize=14)
    axes[1].set_xlabel("Dimension", fontsize=12)
    axes[1].set_ylabel("CPU Time (s)", fontsize=12)
    axes[1].legend(fontsize=12)

    plt.tight_layout()
    plt.savefig("images/performance_Python.png", dpi=150)
    print(f"\nPlot saved to profiling/images/performance_Python.png")

except ImportError:
    print("\nmatplotlib not installed — skipping plot generation.")
    print("Install with: pip install matplotlib")
