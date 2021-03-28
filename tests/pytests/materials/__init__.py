from .TestMaterial import TestMaterial
from .TestElasticity import TestElasticity
from .TestIncompressibleElasticity import TestIncompressibleElasticity
from .TestPoroelasticity import TestPoroelasticity

from .TestRheologies import (
    TestIsotropicLinearElasticity,
    TestIsotropicLinearMaxwell,
    TestIsotropicLinearGenMaxwell,
    TestIsotropicPowerLaw,
    TestIsotropicLinearIncompElasticity,
    TestIsotropicLinearPoroelasticity,
    )

from .TestAuxiliarySubfields import (
    TestAuxSubfieldsElasticity,
    TestAuxSubfieldsIsotropicLinearElasticity,
    TestAuxSubfieldsIsotropicLinearMaxwell,
    TestAuxSubfieldsIsotropicLinearGenMaxwell,
    TestAuxSubfieldsIsotropicPowerLaw,
    TestAuxSubfieldsPoroelasticity,
)

from .TestDerivedSubfields import (
    TestDerivedSubfieldsElasticity,
)

from .TestHomogeneous import TestHomogeneous


def test_classes():
    return [
        TestMaterial,
        TestElasticity,
        TestIncompressibleElasticity,
        TestPoroelasticity,

        TestIsotropicLinearElasticity, 
        TestIsotropicLinearMaxwell,
        TestIsotropicLinearGenMaxwell,
        TestIsotropicPowerLaw,
        TestIsotropicLinearIncompElasticity,
        TestIsotropicLinearPoroelasticity,

        TestAuxSubfieldsElasticity,
        TestAuxSubfieldsIsotropicLinearElasticity,
        TestAuxSubfieldsIsotropicLinearMaxwell,
        TestAuxSubfieldsIsotropicLinearGenMaxwell,
        TestAuxSubfieldsIsotropicPowerLaw,
        TestAuxSubfieldsPoroelasticity,

        TestDerivedSubfieldsElasticity,

        TestHomogeneous,
    ]


# End of file
