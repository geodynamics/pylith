from . import (
    TestDistributor,
    TestField,
    TestMesh,
    TestRefineUniform,
    TestSubfield,
)


def test_modules():
    modules = [
        TestDistributor,
        TestField,
        TestMesh,
        TestRefineUniform,
        TestSubfield,
    ]
    return modules


# End of file
