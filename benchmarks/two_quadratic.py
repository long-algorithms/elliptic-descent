"""Verified type-2+2 relative descents, with independent eclib counts.

Run with sage -python. One JSON line per fixture/repetition. No Magma timing
comparison or full v0.4 implementation is claimed.
"""
import argparse
import json
import time

from foundation import metadata
from export_two_quadratic_arithmetic import FIXTURES
from elliptic_descent import __version__, relative_four_descent

EXPECTED_MODELS = {
    "two_real_fields": 1, "isomorphic_real_fields": 2, "nonmonic_positive": 2,
    "two_imaginary_fields": 1, "mixed_signature": 1, "affine_positive": 2,
    "real_obstruction": 0, "two_adic_obstruction": 0,
}
POINTS = {
    "two_real_fields": (1, 1, 0), "isomorphic_real_fields": (1, 1, 0),
    "nonmonic_positive": (1, 2, 1), "two_imaginary_fields": (1, 1, 0),
    "mixed_signature": (1, 1, 0), "affine_positive": (1, 4, 1),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", choices=["all", *FIXTURES], default="all")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    from sage.all import EllipticCurve, QQ
    from sage.env import SAGE_VERSION
    names = tuple(FIXTURES) if args.fixture == "all" else (args.fixture,)
    for name in names:
        for repetition in range(args.repeats):
            start = time.perf_counter()
            result = relative_four_descent(FIXTURES[name], proof=True)
            solver_seconds = time.perf_counter()-start
            verified, independent_dimension, verification_seconds = False, None, None
            if args.verify:
                start = time.perf_counter()
                verified = result.certified and result.verify_evidence()
                if not verified or len(result.unoriented_coverings) != EXPECTED_MODELS[name]:
                    raise ArithmeticError(f"descent verification failed: {name}")
                if name in POINTS:
                    result.source.evaluate(POINTS[name])
                    E = EllipticCurve([QQ(str(c)) for c in result.source.ainvs])
                    independent_dimension = int(E.mwrank_curve(verbose=False).selmer_rank())
                    if len(result.coverings) != 2**independent_dimension:
                        raise ArithmeticError(f"independent eclib count disagrees: {name}")
                verification_seconds = time.perf_counter()-start
            local = []
            for candidate in result.candidates:
                for report in candidate.local.finite:
                    certificate = report.certificate
                    local.append({"label": candidate.label, "prime": report.place,
                                  "status": report.status.value, "nodes": report.nodes,
                                  "depth": report.depth,
                                  "certificate": certificate.kind if certificate else None,
                                  "model_changes": certificate.change.steps
                                      if certificate and certificate.kind == "model_change" else 0})
            print(json.dumps({
                **metadata(), "package_version": __version__, "sage": SAGE_VERSION,
                "fixture": name, "repetition": repetition, "f_ascending": FIXTURES[name],
                "scope": result.scope, "factorization_type": result.factorization_type,
                "proof": result.arithmetic_proof.value, "assumptions": result.assumptions,
                "support": result.support, "norm_dimension": result.norm_dimension,
                "scalar_image_dimension": len(result.scalar_basis),
                "total_model_classes": result.total_model_classes,
                "processed": len(result.candidates), "unresolved": len(result.unresolved),
                "unoriented_coverings": len(result.unoriented_coverings),
                "oriented_coverings": len(result.coverings),
                "complete": result.complete, "certified": result.certified,
                "verified": verified, "eclib_selmer_dimension": independent_dimension,
                "source_point": POINTS.get(name), "max_nodes": None, "max_depth": None,
                "max_candidates": None, "solver_seconds": solver_seconds,
                "verification_seconds": verification_seconds, **dict(result.timings),
                "local_decisions": local,
            }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
