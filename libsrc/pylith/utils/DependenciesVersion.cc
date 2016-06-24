// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "PetscVersion.hh" // Implementation of class methods

#include "mpi.h"
#include "H5public.h"

// ----------------------------------------------------------------------
#if defined(MPICH_VERSION)
const char* pylith::utils::DependenciesVersion::_mpiVersion = MPICH_VERSION;
const char* pylith::utils::DependenciesVersion::_mpiImplementation = "MPICH";
#else
#if defined(OPENMPI_VERSION)
const char* pylith::utils::DependenciesVersion::_mpiVersion = OPENMPI_VERSION;
const char* pylith::utils::DependenciesVersion::_mpiImplementation = "OpenMPI";
#else
const char* pylith::utils::DependenciesVersion::_mpiVersion = "unknown";
const char* pylith::utils::DependenciesVersion::_mpiImplementation = "unknown";
#endif // OPENMPI
#endif // MPICH

#if defined(NETCDF4_VERSION)
const char* pylith::utils::DependenciesVersion::_netcdfVersion = NETCDF4_VERSION;
#else
const char* pylith::utils::DependenciesVersion::_netcdfVersion = "unknown";
#endif

#if defined(H5_VERS_INFO)
const char* pylith::utils::DependenciesVersion::_hdf5Version = H5_VERS_INFO;
#else
const char* pylith::utils::DependenciesVersion::_hdf5Version = H5_VERS_INFO;
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
