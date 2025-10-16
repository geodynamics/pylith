from . import (
    TestInitializer,
    TestSerial,
    TestParallel,
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
        TestParallel,
        TestConvert,
        TestInitializePhase,
        TestMeshReader,
        TestMeshWriter,
        TestMeshReordering,
        TestMeshRefiner,
        TestMeshDistributor,
        TestMeshInsertInterfaces,
    ]
