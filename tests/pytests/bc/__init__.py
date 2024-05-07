from . import (
    TestDirichletTimeDependent,
    TestNeumannTimeDependent,
    TestAbsorbingDampers,
    TestAuxiliarySubfields,
    TestZeroDB,
    )


def test_modules():
    return [
        TestDirichletTimeDependent,
        TestNeumannTimeDependent,
        TestAbsorbingDampers,
        TestAuxiliarySubfields,
        TestZeroDB,
    ]


# End of file
