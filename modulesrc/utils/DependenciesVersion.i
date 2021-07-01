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
