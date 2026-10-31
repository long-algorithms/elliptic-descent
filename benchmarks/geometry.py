"""Reproducible exact-geometry microbenchmark; not a full descent benchmark."""
import argparse
import json
import platform
import statistics
import time

import sympy as sp

from elliptic_descent import construct_four_cover

FIXTURES = {
    "irreducible": ((1, -1, 0, 0, 1), (0, -1, 0, 0)),
    "nonmonic": ((1, -1, 0, 0, 4), (1, 0, 0, 0)),
    "split": ((4, 0, -5, 0, 1), (1, 0, 0, 0)),
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--repeats", type=int, default=5)
    p.add_argument("--verify", action="store_true")
    args = p.parse_args()
    if args.repeats < 1:
        p.error("--repeats must be positive")
    for name, (f, xi) in FIXTURES.items():
        durations = []
        for _ in range(args.repeats):
            start = time.perf_counter_ns()
            D = construct_four_cover(f, xi)
            if args.verify:
                assert D.check_algebra()
            durations.append((time.perf_counter_ns() - start) / 1e6)
        print(json.dumps({
            "fixture": name, "repeats": args.repeats, "verified": args.verify,
            "median_ms": statistics.median(durations), "samples_ms": durations,
            "python": platform.python_version(), "sympy": sp.__version__,
            "machine": platform.machine(),
        }, sort_keys=True))


if __name__ == "__main__":
    main()
