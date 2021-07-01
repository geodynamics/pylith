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
 * @file modulesrc/utils/PetscVersion.i
 *
 * @brief C++ object for PETSc version information.
 */

namespace pylith {
  namespace utils {

    class PetscVersion
    { // PetscVersion

      // PUBLIC MEMBERS ///////////////////////////////////////////////////////
    public :

      /// Default constructor.
      PetscVersion(void);
      
      /// Default destrictor.
      ~PetscVersion(void);
      
      /** Is source from a release?
       *
       * @returns True if source code comes from a release?
       */
      static
      bool isRelease(void);
      
      /** Get version number.
       *
       * @returns Version number.
       */
      static
      const char* version(void);
      
      /** Get GIT revision.
       *
       * @returns GIT revision.
       */
      static
      const char* gitRevision(void);
      
      /** Get date of GIT revision.
       *
       * @returns Date of GIT revision.
       */
      static
      const char* gitDate(void);
      
      /** Get GIT branch.
       *
       * @returns GIT branch.
       */
      static
      const char* gitBranch(void);
      
      /** Get PETSC_DIR.
       *
       * @returns PETSC_DIR.
       */
      static
      const char* petscDir(void);
      
      /** Get PETSC_ARCH.
       *
       * @returns PETSC_ARCH.
       */
      static
      const char* petscArch(void);
      
    }; // PetscVersion
    
  } // utils
} // pylith

// End of file 
