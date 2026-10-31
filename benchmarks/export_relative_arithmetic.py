"""Export exact arithmetic fixtures for independently replaying local algorithms.

The JSON preserves Sage-computed bases and their provenance. Loading this file
alone is not a new class-group certificate or a complete descent calculation.
"""
import argparse
import json
from pathlib import Path
from time import perf_counter

from elliptic_descent.backends.sage import sage_square_classes
from elliptic_descent.quotient import quotient_fibre

FIXTURES = {
    "x4-x+1": (1, -1, 0, 0, 1),
    "handbook-positive": (172, -120, 20, 12, -7),
    "571a1-first": (-64, 164, -52, -68, -11),
    "571a1-second": (-7, -68, -142, 112, -19),
    "571a1-third": (-3, 52, -232, 60, -4),
    "x4+1": (1, 0, 0, 0, 1),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for name, f in FIXTURES.items():
        print(f"Computing certified arithmetic: {name}", flush=True)
        start = perf_counter()
        context = sage_square_classes(f, proof=True)
        fibre = context.space.norm_fibre()
        scalars = tuple(context.coordinates((c, 0, 0, 0))
                        for c in (-1, *context.space.support))
        quotient, rows = quotient_fibre(fibre, scalars)
        record = {
            "fixture": name,
            "f": [str(c) for c in f],
            "generators": [[str(c) for c in g] for g in context.space.generators],
            "support": list(context.space.support),
            "scalar_vectors": scalars,
            "norm_particular": fibre.particular,
            "norm_basis": fibre.basis,
            "quotient_particular": quotient.particular,
            "quotient_basis": quotient.basis,
            "norm_count": fibre.cardinality,
            "quotient_count": quotient.cardinality,
            "scalar_rank": len(rows),
            "backend": context.space.backend,
            "arithmetic_proof": context.space.arithmetic_proof.value,
            "relative_descent_complete": False,
            "arithmetic_seconds": perf_counter()-start,
        }
        path = args.output / (name + ".json")
        path.write_text(json.dumps(record, indent=2, sort_keys=True)+"\n", encoding="utf-8")
        print(f"{name}: {len(context.space.generators)} generators, "
              f"{fibre.cardinality} norm candidates, {quotient.cardinality} quotient classes; "
              f"{record['arithmetic_seconds']:.3f} seconds", flush=True)


if __name__ == "__main__":
    main()
