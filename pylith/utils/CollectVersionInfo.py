# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/utils/CollectVersionInfo.py
#
# @brief Python CollectVersionInfo object to collect version infofmation for PyLith
# and its dependencies.

from pythia.pyre.components.Component import Component

import pylith.utils.utils as utils

import platform
import sys


class CollectVersionInfo(Component):
    """Python CollectVersionInfo object to collect version information for PyLith
    and its dependencies.
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self):
        """Constructor.
        """
        Component.__init__(self, name="collectversioninfo", facility="collectversioninfo")
        return

    @classmethod
    def asString(cls):
        info = cls._collect()
        s = "Platform:\n" \
            "    Hostname: %(hostname)s\n" \
            "    Operating system: %(os)s\n" \
            "    Kernel: %(kernel)s\n" \
            "    Version: %(version)s\n" \
            "    Machine: %(machine)s\n" \
            "    Processor: %(processor)s\n" \
            % info["platform"]

        version = info["version"]

        # PyLith
        s += "\nPyLith\n"
        if version["pylith"]["isRelease"]:
            s += "    Release v%(version)s\n" % version["pylith"]
        else:
            s += "    Configured on %(gitDate)s, GIT branch: %(gitBranch)s, revision: %(gitRevision)s, hash: %(gitHash)s\n" % version["pylith"]

        # PETSc
        s += "\nPETSc\n"
        if version["petsc"]["isRelease"]:
            s += "    Release v%(version)s\n" % version["petsc"]
        else:
            s += "    Configured on %(gitDate)s, GIT branch: %(gitBranch)s, revision: %(gitRevision)s\n" % version["petsc"]

        # Spatialdata
        s += "\nSpatialdata\n"
        if version["spatialdata"]["isRelease"]:
            s += "    Release v%(version)s\n" % version["spatialdata"]
        else:
            s += "    Configured on %(gitDate)s, GIT branch: %(gitBranch)s, revision: %(gitRevision)s, hash: %(gitHash)s\n" % version["spatialdata"]

        # MPI
        s += "\nMPI standard: %(standard)s, implementation: %(implementation)s, version: %(version)s\n" % version["mpi"]

        # HDF5
        s += "HDF5 version: %(version)s\n" % version["hdf5"]

        # NetCDF
        s += "NetCDF version: %(version)s\n" % version["netcdf"]

        # Proj
        s += "Proj version: %(version)s\n" % version["proj"]

        # Python
        s += "\nPython\n" \
            "    v%(version)s of %(implementation)s compiled with %(compiler)s\n" % version["python"]
        for (pname, pver) in version["python"]["modules"].items():
            s += "    %s: v%s from %s\n" % (pname, pver["version"], pver["location"])
        return s

    @classmethod
    def asDict(cls):
        info = cls._collect()
        return info

    # PRIVATE METHODS ////////////////////////////////////////////////////

    @classmethod
    def _collect(cls):
        """Collect version infoformation.
        """
        info = {
            "platform": cls._collectPlatform(),
            "version": cls._collectVersion(),
        }
        return info

    @classmethod
    def _collectPlatform(cls):
        (os, hostname, kernel, version, machine, processor) = platform.uname()
        info = {
            "hostname": hostname,
            "os": os,
            "kernel": kernel,
            "version": version,
            "machine": machine,
            "processor": processor,
        }
        return info

    @classmethod
    def _collectVersion(cls):
        info = {
            "pylith": cls._collectVersionPyLith(),
            "python": cls._collectVersionPython(),
            "petsc": cls._collectVersionPetsc(),
            "mpi": cls._collectVersionMPI(),
            "hdf5": cls._collectVersionHDF5(),
            "netcdf": cls._collectVersionNetCDF(),
            "spatialdata": cls._collectVersionSpatialdata(),
            "proj": cls._collectVersionProj(),
        }
        return info

    @staticmethod
    def _collectVersionPyLith():
        v = utils.PylithVersion()
        if v.isRelease():
            info = {
                "isRelease": True,
                "version": v.version()
            }
        else:
            info = {
                "isRelease": False,
                "gitDate": v.gitDate(),
                "gitBranch": v.gitBranch(),
                "gitRevision": v.gitRevision(),
                "gitHash": v.gitHash(),
            }
        return info

    @staticmethod
    def _collectVersionSpatialdata():
        import spatialdata.utils.utils as utils
        v = utils.SpatialdataVersion()
        if v.isRelease():
            info = {
                "isRelease": True,
                "version": v.version()
            }
        else:
            info = {
                "isRelease": False,
                "gitDate": v.gitDate(),
                "gitBranch": v.gitBranch(),
                "gitRevision": v.gitRevision(),
                "gitHash": v.gitHash(),
            }
        return info

    @staticmethod
    def _collectVersionPetsc():
        v = utils.PetscVersion()
        if v.isRelease():
            info = {
                "isRelease": True,
                "version": v.version()
            }
        else:
            info = {
                "isRelease": False,
                "gitDate": v.gitDate(),
                "gitBranch": v.gitBranch(),
                "gitRevision": v.gitRevision(),
            }
        info["petscDir"] = v.petscDir()
        info["petscArch"] = v.petscArch()
        return info

    @classmethod
    def _collectVersionPython(cls):
        info = {
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
            "compiler": platform.python_compiler(),
            "modules": {},
        }
        pkgs = ("numpy", "spatialdata", "h5py", "netCDF4", "pythia")
        for pkg in pkgs:
            ver, loc = cls._getPackageVersion(pkg)
            info["modules"][pkg] = {
                "version": ver,
                "location": loc,
            }
        return info

    @staticmethod
    def _collectVersionMPI():
        v = utils.DependenciesVersion()
        info = {
            "standard": v.mpiStandard(),
            "implementation": v.mpiImplementation(),
            "version": v.mpiVersion(),
        }
        return info

    @staticmethod
    def _collectVersionHDF5():
        v = utils.DependenciesVersion()
        info = {
            "version": v.hdf5Version(),
        }
        return info

    @staticmethod
    def _collectVersionNetCDF():
        v = utils.DependenciesVersion()
        info = {
            "version": v.netcdfVersion(),
        }
        return info

    @staticmethod
    def _collectVersionProj():
        import spatialdata.utils.utils as utils
        v = utils.SpatialdataVersion()
        info = {
            "version": v.projVersion(),
        }
        return info

    @staticmethod
    def _getPackageVersion(name):
        import os
        m = None
        location = None
        version = None
        try:
            m = __import__(name)
            location = os.path.split(m.__file__)[0]
            version = m.__version__
        except ImportError:
            version = "not found"
            location = "--"
        except AttributeError:
            if version is None:
                version = "unknown"
            if location is None:
                location = "unknown"
        return (version, location)


# End of file
