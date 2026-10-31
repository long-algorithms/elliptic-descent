"""Sage checks of exact split-torsion presentation selection."""
import os
import pytest

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE') == '1':
    import sage.all as sage
else:
    sage=pytest.importorskip('sage.all',reason='requires SageMath runtime')

from elliptic_descent import two_selmer


def test_split_torsion_presentation_selection_preserves_labels_and_all_support():
    from elliptic_descent.fibre_presentations import fibre_representative
    from elliptic_descent.arithmetic import quartic_support
    S=two_selmer((0,0,0,-289,0))
    torsion=S.torsion_data()
    for entry in torsion.entries:
        C=fibre_representative(S,entry.coordinates,torsion)
        assert C.certified and C.verify()
        assert C.ambient_label==S.ambient_label(entry.coordinates)
        if entry.coordinates:
            assert quartic_support(C.quartic.f)==(2,17)
    # The selected models must still pass actual relative descent and all
    # equivalence checks; the separate eight-curve benchmark exercises that.
