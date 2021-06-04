from .TestDistributor import TestDistributor
from .TestField import TestField
from .TestMesh import TestMesh
from .TestMeshGenerator import TestMeshGenerator
from .TestMeshImporter import TestMeshImporter
from .TestMeshRefiner import TestMeshRefiner
from .TestRefineUniform import TestRefineUniform
from .TestReverseCuthillMcKee import TestReverseCuthillMcKee
from .TestSubfield import TestSubfield


def test_classes():
    classes = [
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
    return classes


# End of file
