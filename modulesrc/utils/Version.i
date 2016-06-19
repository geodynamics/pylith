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
 * @file modulesrc/utils/Version.hh
 *
 * @brief C++ object for PyLith version information.
 */

namespace pylith {
  namespace utils {

    class pylith::utils::Version
    { // Version

      // PUBLIC MEMBERS ///////////////////////////////////////////////////////
    public :

      /// Default constructor.
      Version(void);
      
      /// Default destrictor.
      ~Version(void);

      /** Get Pylith version number.
       *
       * @returns PyLith version numer.
       */
      static 
	const char* version(void);

    }; // Version

  } // utils
} // pylith

// End of file 
