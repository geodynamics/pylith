// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================

/**
 * @file modulesrc/utils/PylithVersion.i
 *
 * @brief C++ object for PyLith version information.
 */

namespace pylith {
  namespace utils {

    class DependenciesVersion
    { // DependenciesPylithVersion

      // PUBLIC MEMBERS ///////////////////////////////////////////////////////
    public :

      /// Default constructor.
      DependenciesVersion(void);

      /// Default destrictor.
      ~DependenciesVersion(void);

      /** Get MPI version number.
       *
       * @returns MPI version number.
       */
      static
      const char* mpiVersion(void);
      
      /** Get MPI implemenation info (OpenMPI, MPICH, etc).
       *
       * @returns MPI implementation info.
       */
      static
      const char* mpiImplementation(void);
      
      /** Get MPI standard version info.
       *
       * @returns MPI standard version info.
       */
      static
      const char* mpiStandard(void);

      /** Get NetCDF version number.
       *
       * @returns NetCDF version number.
       */
      static
      const char* netcdfVersion(void);
      
      /** Get HDF5 version number.
       *
       * @returns HDF5 version number.
       */
      static
      const char* hdf5Version(void);      
      
    }; // DependenciesPylithVersion
    
  } // utils
} // pylith

// End of file 
