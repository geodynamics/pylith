from . import (
    TestCollectVersionInfo,
    TestConstants,
    TestEmptyBin,
    TestNullComponent,
    TestDumpParameters,
    TestDumpParametersAscii,
    TestDumpParametersJson,
    TestDependenciesVersion,
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
        TestDependenciesVersion,
        TestPylithVersion,
        TestProfiling,
    ]
    return modules


# End of file
