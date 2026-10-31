"""CRT and factorwise square-class benchmarks; NOT a complete descent benchmark."""
import argparse
import json
import time

from foundation import metadata
from elliptic_descent import EtaleQuarticAlgebra, __version__
from elliptic_descent.geometry import _norm

FIXTURES = {
    "irreducible": (1, -1, 0, 0, 1),
    "cubic_linear": (-3, 4, -1, -3, 1),
    "two_quadratics": (6, 0, -5, 0, 1),
    "quadratic_two_linear": (2, 0, -3, 0, 1),
    "split": (4, 0, -5, 0, 1),
    "isomorphic_quadratics": (16, 0, -10, 0, 1),
    "class_group_two_torsion": (-5, 0, 4, 0, 1),
    "nonmonic": (12, 0, -10, 0, 2),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--sage", action="store_true")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    for name, f in FIXTURES.items():
        for repetition in range(args.repeats):
            start = time.perf_counter()
            algebra = EtaleQuarticAlgebra(f)
            # Force the lazy factorization and both CRT maps into this stage.
            parts = algebra.to_components((1, 2, 3, 4))
            roundtrip = algebra.from_components(parts)
            construction_seconds = time.perf_counter()-start
            context = None
            arithmetic_seconds = None
            if args.sage:
                from elliptic_descent.backends.sage_etale import sage_etale_square_classes
                start = time.perf_counter()
                context = sage_etale_square_classes(algebra, proof=True)
                arithmetic_seconds = time.perf_counter()-start
            verification_seconds = None
            if args.verify:
                start = time.perf_counter()
                if (roundtrip != (1, 2, 3, 4) or not algebra.verify_decomposition()
                        or algebra.norm(roundtrip) != _norm(roundtrip, algebra.f)):
                    raise ArithmeticError(f"CRT verification failed: {name}")
                if context is not None:
                    if context.norm_fibre() != context.space.norm_fibre():
                        raise ArithmeticError(f"norm-fibre verification failed: {name}")
                    for j, g in enumerate(context.space.generators):
                        if (context.coordinates(g) != 1 << j
                                or context.coordinates(context.representative(1 << j)) != 1 << j):
                            raise ArithmeticError(f"basis verification failed: {name}")
                verification_seconds = time.perf_counter()-start
            record = {
                **metadata(), "package_version": __version__, "fixture": name,
                "repetition": repetition, "f_ascending": f,
                "factor_degrees": algebra.degrees, "verified": args.verify,
                "scope": "quartic CRT arithmetic only" if context is None else "A(S,2) and norm fibre only",
                "relative_descent_complete": False,
                "construction_seconds": construction_seconds,
                "arithmetic_seconds": arithmetic_seconds,
                "verification_seconds": verification_seconds,
            }
            if context is not None:
                record.update({
                    "backend": context.space.backend, "support": context.space.support,
                    "factor_dimensions": context.factor_dimensions,
                    "square_class_dimension": context.dimension,
                    "norm_fibre_cardinality": context.norm_fibre().cardinality,
                    "proof": context.space.arithmetic_proof.value,
                    "assumptions": context.space.assumptions,
                })
            print(json.dumps(record, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
