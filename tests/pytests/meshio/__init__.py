from .TestDataWriter import TestDataWriter
from .TestDataWriterVTK import TestDataWriterVTK
from .TestMeshIOAscii import TestMeshIOAscii
from .TestOutputObserver import TestOutputObserver
from .TestOutputPhysics import TestOutputPhysics
from .TestOutputSoln import TestOutputSoln
from .TestOutputSolnDomain import TestOutputSolnDomain
from .TestOutputSolnBoundary import TestOutputSolnBoundary
from .TestOutputSolnPoints import TestOutputSolnPoints
from .TestOutputTrigger import TestOutputTrigger
from .TestOutputTriggerStep import TestOutputTriggerStep
from .TestOutputTriggerTime import TestOutputTriggerTime
from .TestPointsList import TestPointsList


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


def test_classes():
    classes = [
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
        from .TestMeshIOCubit import TestMeshIOCubit
        classes += [
            TestMeshIOCubit,
        ]
    if has_h5py():
        from .TestDataWriterHDF5 import TestDataWriterHDF5
        from .TestDataWriterHDF5Ext import TestDataWriterHDF5Ext
        from .TestXdmf import TestXdmf
        classes += [
            TestDataWriterHDF5,
            TestDataWriterHDF5Ext,
            TestXdmf,
        ]
    return classes


# End of file
