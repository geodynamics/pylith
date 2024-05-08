from . import (
    TestCollectVersionInfo,
    TestConstants,
    TestEmptyBin,
    TestNullComponent,
    TestDumpParameters,
    TestDumpParametersAscii,
    TestDumpParametersJson,
    TestEventLogger,
    TestPetscManager,
    TestDependenciesVersion,
    TestPetscVersion,
    TestPylithVersion,
    TestProfiling,
)


def test_modules():
    modules = [
        TestCollectVersionInfo,
        TestConstants,
        TestEmptyBin,
        TestNullComponent,
        TestDumpParameters,
        TestDumpParametersAscii,
        TestDumpParametersJson,
        TestEventLogger,
        TestPetscManager,
        TestDependenciesVersion,
        TestPetscVersion,
        TestPylithVersion,
        TestProfiling,
    ]
    return modules


# End of file
