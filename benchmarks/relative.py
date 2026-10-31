"""Reproducible relative-descent benchmarks; run with sage -python.

Each JSON line records one fixture and repetition. Verification is separate
from the solver timing. These are reference implementations, not Magma timing
comparisons; the published counts use the explicitly un-oriented quotient.
"""
import argparse
import json
import platform
import time

import sympy as sp

from foundation import metadata
from export_relative_arithmetic import FIXTURES
from elliptic_descent import __version__
from elliptic_descent.relative import relative_four_descent

EXPECTED_MODELS = {"x4-x+1": 2, "handbook-positive": 1,
                   "571a1-first": 0, "571a1-second": 0, "571a1-third": 0,
                   "x4+1": 2}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", choices=["all", *FIXTURES], default="all")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--max-nodes", type=int)
    parser.add_argument("--max-depth", type=int)
    parser.add_argument("--max-candidates", type=int)
    args = parser.parse_args()
    if args.repeats < 1 or any(v is not None and v < 0 for v in
                             (args.max_nodes,args.max_depth,args.max_candidates)):
        parser.error("repeats must be positive; limits must be nonnegative")
    from sage.env import SAGE_VERSION
    names = tuple(FIXTURES) if args.fixture == "all" else (args.fixture,)
    for name in names:
        for repetition in range(args.repeats):
            start = time.perf_counter()
            result = relative_four_descent(FIXTURES[name], proof=True,
                         max_candidates=args.max_candidates, max_depth=args.max_depth,
                         max_nodes=args.max_nodes)
            solver_seconds = time.perf_counter()-start
            verification_seconds, verified = None, False
            if args.verify:
                start = time.perf_counter()
                verified = result.certified and result.verify_evidence()
                verification_seconds = time.perf_counter()-start
                if not verified or len(result.unoriented_coverings) != EXPECTED_MODELS[name]:
                    raise RuntimeError(f"incomplete or incorrect benchmark result: {name}")
            decisions = []
            coefficient_bits = 0
            for candidate in result.candidates:
                for q in candidate.model.quadrics:
                    for c in sp.Poly(q).coeffs():
                        coefficient_bits = max(coefficient_bits, abs(int(c.p)).bit_length(),
                                                int(c.q).bit_length())
                for r in candidate.local.finite:
                    c = r.certificate
                    decisions.append({"label": candidate.label, "prime": r.place,
                                      "status": r.status.value, "nodes": r.nodes,
                                      "depth": r.depth, "certificate": c.kind if c else None,
                                      "model_changes": c.change.steps if c and c.kind == "model_change" else 0})
            print(json.dumps({**metadata(), "fixture": name, "repetition": repetition,
                              "package_version": __version__, "sage": SAGE_VERSION,
                              "machine": platform.machine(), "f_ascending": FIXTURES[name],
                              "scope": result.scope, "proof": result.arithmetic_proof.value,
                              "assumptions": result.assumptions, "support": result.support,
                              "norm_dimension": result.norm_dimension,
                              "scalar_image_dimension": len(result.scalar_basis),
                              "total_model_classes": result.total_model_classes,
                              "processed": len(result.candidates), "unresolved": len(result.unresolved),
                              "unoriented_coverings": len(result.unoriented_coverings),
                              "oriented_coverings": len(result.coverings),
                              "complete": result.complete, "certified": result.certified,
                              "verified": verified, "max_coefficient_bits": coefficient_bits,
                              "max_nodes": args.max_nodes, "max_depth": args.max_depth,
                              "max_candidates": args.max_candidates,
                              "solver_seconds": solver_seconds, "verification_seconds": verification_seconds,
                              **dict(result.timings), "local_decisions": decisions}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
