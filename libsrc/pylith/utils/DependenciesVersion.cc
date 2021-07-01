// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "DependenciesVersion.hh" // Implementation of class methods

#include "mpi.h"
#include "H5pubconf.h"

// ----------------------------------------------------------------------
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)


#if defined(MPICH_VERSION)
const char* pylith::utils::DependenciesVersion::_mpiVersion = MPICH_VERSION;
const char* pylith::utils::DependenciesVersion::_mpiImplementation = "MPICH";
#else
#if defined(OMPI_MAJOR_VERSION)
#define PYLITH_OPENMPI_VERSION STR(OMPI_MAJOR_VERSION) "." STR(OMPI_MINOR_VERSION) "." STR(OMPI_RELEASE_VERSION)

const char* pylith::utils::DependenciesVersion::_mpiVersion = PYLITH_OPENMPI_VERSION;
const char* pylith::utils::DependenciesVersion::_mpiImplementation = "OpenMPI";
#else
const char* pylith::utils::DependenciesVersion::_mpiVersion = "unknown";
const char* pylith::utils::DependenciesVersion::_mpiImplementation = "unknown";
#endif // OPENMPI
#endif // MPICH
#define PYLITH_MPI_STANDARD STR(MPI_VERSION) "." STR(MPI_SUBVERSION)
const char* pylith::utils::DependenciesVersion::_mpiStandard = PYLITH_MPI_STANDARD;

#if defined(NETCDF4_VERSION)
const char* pylith::utils::DependenciesVersion::_netcdfVersion = NETCDF4_VERSION;
#else
const char* pylith::utils::DependenciesVersion::_netcdfVersion = "unknown";
#endif

#if defined(H5_VERSION)
const char* pylith::utils::DependenciesVersion::_hdf5Version = H5_VERSION;
#else
const char* pylith::utils::DependenciesVersion::_hdf5Version = H5_VERSION;
#endif

// ----------------------------------------------------------------------
// Default constructor.
pylith::utils::DependenciesVersion::DependenciesVersion(void)
{}

// ----------------------------------------------------------------------
// Default destrictor.
pylith::utils::DependenciesVersion::~DependenciesVersion(void)
{}

// ----------------------------------------------------------------------
// Get MPI version number.
const char*
pylith::utils::DependenciesVersion::mpiVersion(void)
{ // mpiVersion
  return _mpiVersion;
} // mpiVersion

// ----------------------------------------------------------------------
// Get MPI version number.
const char*
pylith::utils::DependenciesVersion::mpiImplementation(void)
{ // mpiImplementation
  return _mpiImplementation;
} // mpiImplementation

// ----------------------------------------------------------------------
// Get MPI standard version number.
const char*
pylith::utils::DependenciesVersion::mpiStandard(void)
{ // mpiStandard
  return _mpiStandard;
} // mpiStandard

// ----------------------------------------------------------------------
// Get NetCDF version number.
const char*
pylith::utils::DependenciesVersion::netcdfVersion(void)
{ // mpiVersion
  return _netcdfVersion;
} // netcdfVersion

// ----------------------------------------------------------------------
// Get HDF5 version number.
const char*
pylith::utils::DependenciesVersion::hdf5Version(void)
{ // hdf5Version
  return _hdf5Version;
} // hdf5Version

  
// End of file 
