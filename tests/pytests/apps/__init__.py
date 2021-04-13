from .TestPetscApplication import TestPetscApplication
from .TestPyLithApp import TestPyLithApp
from .TestEqInfoApp import TestEqInfoApp


def test_classes():
    return [
        TestPetscApplication,
        TestPyLithApp,
        TestEqInfoApp,
    ]


# End of file
