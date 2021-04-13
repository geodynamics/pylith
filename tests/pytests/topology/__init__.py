from .TestDistributor import TestDistributor
from .TestField import TestField
from .TestFields import TestFields
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
        TestFields,
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
