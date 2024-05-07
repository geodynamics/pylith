from . import (
    TestDistributor,
    TestField,
    TestMesh,
    TestMeshGenerator,
    TestMeshImporter,
    TestMeshRefiner,
    TestRefineUniform,
    TestReverseCuthillMcKee,
    TestSubfield,
)


def test_modules():
    modules = [
        TestDistributor,
        TestField,
        TestMesh,
        TestMeshGenerator,
        TestMeshImporter,
        TestMeshRefiner,
        TestRefineUniform,
        TestReverseCuthillMcKee,
        TestSubfield,
    ]
    return modules


# End of file
