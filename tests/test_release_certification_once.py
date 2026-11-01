"""Independent-audit certification contract; not a mathematical stand-in."""
import importlib.util
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_independent_audit_certifies_once_per_invocation_without_memoization():
    # This is an orchestration regression, not a replacement for Sage arithmetic.
    sys.path.insert(0, str(ROOT/'benchmarks'))
    try:
        spec = importlib.util.spec_from_file_location('independent_mw_once', ROOT/'benchmarks/independent_mw.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)

    class Result:
        invariant_factors = (2, 4)
        valid = True
        checks = 0

        def require_group(self, *, certified=False):
            assert certified is False  # checked_group performs the one check below.
            return self

        @property
        def certified(self):
            self.checks += 1
            return self.valid

    result = Result()
    assert module.checked_group(result, (2, 4)) == (result, True)
    assert result.checks == 1
    assert module.checked_group(result, (2, 4)) == (result, True)
    assert result.checks == 2
    result.valid = False
    with pytest.raises(ArithmeticError):
        module.checked_group(result, (2, 4))
    assert result.checks == 3
    result.valid = True
    with pytest.raises(ArithmeticError):
        module.checked_group(result, (4, 4))
    assert result.checks == 4
