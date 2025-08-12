from . import (
    TestInitializer,
    TestSerial,
    TestConvert,
    TestInitializePhase,
    TestMeshReader,
    TestMeshWriter,
    TestMeshReordering,
    TestMeshRefiner,
    TestMeshDistributor,
    TestMeshInsertInterfaces,
)

def test_modules():
    return [
        TestInitializer,
        TestSerial,
        TestConvert,
        TestInitializePhase,
        TestMeshReader,
        TestMeshWriter,
        TestMeshReordering,
        TestMeshRefiner,
        TestMeshDistributor,
        TestMeshInsertInterfaces,
    ]
