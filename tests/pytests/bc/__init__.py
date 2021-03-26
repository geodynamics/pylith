from .TestDirichletTimeDependent import TestDirichletTimeDependent
from .TestNeumannTimeDependent import TestNeumannTimeDependent
from .TestAbsorbingDampers import TestAbsorbingDampers
from .TestAuxiliarySubfields import (
    TestAuxSubfieldsTimeDependent, TestAuxSubfieldsAbsorbingDampers)
from .TestZeroDB import TestZeroDB


def test_classes():
    return [
        TestDirichletTimeDependent,
        TestNeumannTimeDependent,
        TestAuxSubfieldsTimeDependent,
        TestAbsorbingDampers,
        TestAuxSubfieldsAbsorbingDampers,
        TestZeroDB,
    ]


# End of file
