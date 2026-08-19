from . import (
    TestApplication,
    TestEventLogger,
    TestManager,
    TestPetscVersion,
)


def test_modules():
    modules = [
        TestApplication,
        TestEventLogger,
        TestManager,
        TestPetscVersion,
    ]
    return modules


# End of file
