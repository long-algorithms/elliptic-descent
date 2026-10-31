"""Export certified arithmetic for replaying two-quadratic local descents.

The exported bases retain provenance but are not standalone class-group
certificates. No completeness claim follows from loading the JSON alone.
"""
import argparse
import json
from pathlib import Path
from time import perf_counter

from elliptic_descent.backends.sage_etale import sage_etale_square_classes
from elliptic_descent.quotient import quotient_fibre

FIXTURES = {
    "two_real_fields": (6, 0, -5, 0, 1),
    "isomorphic_real_fields": (16, 0, -10, 0, 1),
    "nonmonic_positive": (12, 0, -10, 0, 2),
    "two_imaginary_fields": (2, 0, 3, 0, 1),
    "mixed_signature": (-6, 0, 1, 0, 1),
    "real_obstruction": (-2, 0, -3, 0, -1),
    "affine_positive": (6, 0, 8, 0, 2),
    "two_adic_obstruction": (6, 0, 9, 0, 3),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for name, f in FIXTURES.items():
        start = perf_counter()
        context = sage_etale_square_classes(f, proof=True)
        assert context.algebra.factorization_type == (2, 2)
        fibre = context.norm_fibre()
        scalars = tuple(context.coordinates(context.algebra.scalar(c))
                        for c in (-1, *context.space.support))
        quotient, rows = quotient_fibre(fibre, scalars)
        record = {
            "fixture": name, "f": [str(c) for c in f],
            "generators": [[str(c) for c in g] for g in context.space.generators],
            "support": list(context.space.support), "scalar_vectors": scalars,
            "norm_particular": fibre.particular, "norm_basis": fibre.basis,
            "quotient_particular": quotient.particular, "quotient_basis": quotient.basis,
            "norm_count": fibre.cardinality, "quotient_count": quotient.cardinality,
            "scalar_rank": len(rows), "factor_dimensions": context.factor_dimensions,
            "backend": context.space.backend, "arithmetic_proof": context.space.arithmetic_proof.value,
            "relative_descent_complete": False, "arithmetic_seconds": perf_counter()-start,
        }
        (args.output/(name+".json")).write_text(json.dumps(record, indent=2, sort_keys=True)+"\n")
        print(name, "quotient classes:", quotient.cardinality, flush=True)


if __name__ == "__main__":
    main()
