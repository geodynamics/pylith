from .TestFaultCohesive import TestFaultCohesive
from .TestFaultCohesiveKin import TestFaultCohesiveKin
from .TestKinSrc import (
    TestKinSrc, 
    TestKinSrcConstRate,
    TestKinSrcStep,
    TestKinSrcRamp,
    TestKinSrcBrune,
    TestKinSrcLiuCos,
    TestKinSrcTimeHistory)
from .TestSingleRupture import TestSingleRupture


def test_classes():
    return [
        TestFaultCohesive,
        TestFaultCohesiveKin,
        TestKinSrc,
        TestKinSrcConstRate,
        TestKinSrcStep,
        TestKinSrcRamp,
        TestKinSrcBrune,
        TestKinSrcLiuCos,
        TestKinSrcTimeHistory,
        TestSingleRupture,
    ]


# End of file
