"""JSONL benchmarks for v0.2; timings never stand in for proof status."""
import argparse
import json
import os
import platform
import statistics
import subprocess
import time

import sympy as sp

from elliptic_descent.geometry import U
from elliptic_descent.local import padic_solubility
from elliptic_descent.status import LocalStatus


def metadata():
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL,
                                           text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        revision = os.environ.get("GITHUB_SHA", "unavailable")
    try:
        import resource
        memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    except ImportError:
        memory = None
    return {"revision": revision, "python": platform.python_version(), "sympy": sp.__version__,
            "platform": platform.platform(), "cpu": platform.processor(),
            "process_peak_rss_native_units": memory,
            "rss_units": "bytes on macOS; KiB on Linux; cumulative process peak",
            "cache_mode": "in-process repetitions; not cold-process measurements"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--sage", action="store_true")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    u, v, w, z = U
    fixtures = {
        "obstruction_at_2": ((u*u+v*v+w*w+z*z, u*u+2*v*v+3*w*w+4*z*z), 2),
        "singular_hensel_at_2": ((u*u-17*v*v, w*w-17*z*z), 2),
        "hensel_at_5": ((u*u+v*v+w*w+z*z, u*u+2*v*v+3*w*w+4*z*z), 5),
    }
    for name, (quadrics, p) in fixtures.items():
        samples, verification = [], []
        for _ in range(args.repeats):
            start = time.perf_counter_ns()
            result = padic_solubility(quadrics, p, max_nodes=100_000, max_depth=6)
            samples.append((time.perf_counter_ns()-start)/1e6)
            if args.verify:
                start = time.perf_counter_ns()
                if result.status == LocalStatus.UNKNOWN or not result.verify(quadrics):
                    raise RuntimeError(f"benchmark certificate failed: {name}")
                verification.append((time.perf_counter_ns()-start)/1e6)
        print(json.dumps({**metadata(), "fixture": name, "stage": "local", "prime": p,
                          "status": result.status.value, "certificate": result.certificate.kind,
                          "verified": args.verify, "nodes": result.nodes, "depth": result.depth,
                          "max_nodes": 100_000, "max_depth": 6, "samples_ms": samples,
                          "median_ms": statistics.median(samples),
                          "verification_samples_ms": verification}, sort_keys=True))
    if args.sage:
        from elliptic_descent.backends.sage import sage_square_classes
        from elliptic_descent.foundation import screen_norm_candidates
        start = time.perf_counter_ns()
        context = sage_square_classes((1, -1, 0, 0, 1), proof=True)
        arithmetic_ms = (time.perf_counter_ns()-start)/1e6
        start = time.perf_counter_ns()
        result = screen_norm_candidates(context.space, limit=128, max_nodes=5000,
                                        verify_geometry=args.verify)
        screening_ms = (time.perf_counter_ns()-start)/1e6
        if args.verify:
            for candidate in result.candidates:
                for report in candidate.local.reports:
                    if report.certificate and not report.verify(candidate.covering):
                        raise RuntimeError("Sage pipeline certificate verification failed")
        print(json.dumps({**metadata(), "fixture": "x4-x+1", "stage": "arithmetic_and_screen",
                          "backend": context.space.backend, "proof": "certified",
                          "assumptions": result.arithmetic.assumptions,
                          "support": context.space.support,
                          "basis_dimension": len(context.space.generators),
                          "norm_candidates": result.arithmetic.total_count,
                          "processed": len(result.candidates), "unresolved": len(result.unresolved),
                          "norm_fibre_exhausted": result.arithmetic.norm_fibre_exhausted,
                          "relative_descent_complete": result.relative_descent_complete,
                          "arithmetic_ms": arithmetic_ms, "screening_ms": screening_ms,
                          "geometry_verified": args.verify}, sort_keys=True))


if __name__ == "__main__":
    main()
