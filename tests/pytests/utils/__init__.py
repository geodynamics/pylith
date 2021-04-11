from .TestCollectVersionInfo import TestCollectVersionInfo
from .TestConstants import TestConstants
from .TestEmptyBin import TestEmptyBin
from .TestNullComponent import TestNullComponent
from .TestDumpParameters import TestDumpParameters
from .TestDumpParametersAscii import TestDumpParametersAscii
from .TestDumpParametersJson import TestDumpParametersJson
from .TestEventLogger import TestEventLogger
from .TestPetscManager import TestPetscManager
from .TestDependenciesVersion import TestDependenciesVersion
from .TestPetscVersion import TestPetscVersion
from .TestPylithVersion import TestPylithVersion
from .TestProfiling import TestProfiling


def test_classes():
    classes = [
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
    return classes


# End of file
