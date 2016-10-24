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
 * @file modulesrc/utils/PylithVersion.i
 *
 * @brief C++ object for PyLith version information.
 */

namespace pylith {
  namespace utils {

    class PylithVersion
    { // PylithVersion

      // PUBLIC MEMBERS ///////////////////////////////////////////////////////
    public :

      /// Default constructor.
      PylithVersion(void);
      
      /// Default destrictor.
      ~PylithVersion(void);
      
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
      
      /** Get GIT hash.
       *
       * @returns GIT hash.
       */
      static
      const char* gitHash(void) const;
      
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
      
      
    }; // PylithVersion
    
  } // utils
} // pylith

// End of file 
