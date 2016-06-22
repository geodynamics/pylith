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

/**
 * @file modulesrc/utils/PetscVersion.i
 *
 * @brief C++ object for PETSc version information.
 */

namespace pylith {
  namespace utils {

    class pylith::utils::PetscVersion
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
      bool isRelease(void) const;
      
      /** Get version number.
       *
       * @returns Version number.
       */
      static
      const char* version(void) const;
      
      /** Get GIT revision.
       *
       * @returns GIT revision.
       */
      static
      const char* gitRevision(void) const;
      
      /** Get date of GIT revision.
       *
       * @returns Date of GIT revision.
       */
      static
      const char* gitDate(void) const;
      
      /** Get GIT branch.
       *
       * @returns GIT branch.
       */
      static
      const char* gitBranch(void) const;
      
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
