from . import (
    TestMeshIOAscii,
    TestMeshIOPetsc,
    TestDataWriter,
    TestDataWriterVTK,
    TestOutputObserver,
    TestOutputPhysics,
    TestOutputSoln,
    TestOutputSolnDomain,
    TestOutputSolnBoundary,
    TestOutputSolnPoints,
    TestOutputTrigger,
    TestOutputTriggerStep,
    TestOutputTriggerTime,
    TestPointsList,
)


def has_h5py():
    flag = True
    try:
        import h5py
    except:
        flag = False
    return flag


def has_netcdf():
    flag = True
    try:
        import netCDF4
    except:
        flag = False
    return flag


def test_modules():
    modules = [
        TestMeshIOAscii,
        TestMeshIOPetsc,
        TestDataWriter,
        TestDataWriterVTK,
        TestOutputObserver,
        TestOutputPhysics,
        TestOutputSoln,
        TestOutputSolnDomain,
        TestOutputSolnBoundary,
        TestOutputSolnPoints,
        TestOutputTrigger,
        TestOutputTriggerStep,
        TestOutputTriggerTime,
        TestPointsList,
    ]
    if has_netcdf():
        from . import TestMeshIOCubit
        modules += [
            TestMeshIOCubit,
        ]
    if has_h5py():
        from . import (
            TestDataWriterHDF5,
            TestDataWriterHDF5Ext,
            TestXdmf,
        )
        modules += [
            TestDataWriterHDF5,
            TestDataWriterHDF5Ext,
            TestXdmf,
        ]
    return modules


# End of file
