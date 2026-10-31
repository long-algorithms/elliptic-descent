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
        verify_durations = []
        for _ in range(args.repeats):
            start = time.perf_counter_ns()
            D = construct_four_cover(f, xi)
            durations.append((time.perf_counter_ns() - start) / 1e6)
            if args.verify:
                verification_start = time.perf_counter_ns()
                assert D.check_algebra()
                assert D.check_map_identity()
                assert D.check_smoothness()
                assert D.check_map_basepoint_free()
                verify_durations.append((time.perf_counter_ns() - verification_start) / 1e6)
        print(json.dumps({
            "fixture": name, "repeats": args.repeats, "verified": args.verify,
            "verification": "algebra,map_identity,smoothness,map_basepoint_free" if args.verify else "none",
            "median_ms": statistics.median(durations), "samples_ms": durations,
            "verification_median_ms": (statistics.median(verify_durations) if args.verify else None),
            "verification_samples_ms": verify_durations,
            "python": platform.python_version(), "sympy": sp.__version__,
            "machine": platform.machine(),
        }, sort_keys=True))


if __name__ == "__main__":
    main()
