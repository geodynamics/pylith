from . import (
    TestFaultCohesive,
    TestFaultCohesiveKin,
    TestFaultCohesiveImpulses,
    TestKinSrc, 
    TestSingleRupture,
 )


def test_modules():
    return [
        TestFaultCohesive,
        TestFaultCohesiveKin,
        TestFaultCohesiveImpulses,
        TestKinSrc,
        TestSingleRupture,
    ]


# End of file
