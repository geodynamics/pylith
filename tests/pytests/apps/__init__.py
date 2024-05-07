from . import (TestPetscApplication, TestPyLithApp, TestEqInfoApp)


def test_modules():
    return [
        TestPetscApplication,
        TestPyLithApp,
        TestEqInfoApp,
    ]


# End of file
