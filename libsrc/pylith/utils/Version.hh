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
 * @file libsrc/utils/Version.hh
 *
 * @brief C++ object for PyLith version information.
 */

#if !defined(pylith_utils_version_hh)
#define pylith_utils_version_hh

// Include directives ---------------------------------------------------
#include "utilsfwd.hh" // forward declarations

// Version ----------------------------------------------------------
/** @brief C++ object for getting version info.
 */
class pylith::utils::Version
{ // Version
  friend class TestVersion; // unit testing

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

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  Version(const Version&); ///< Not implemented
  const Version& operator=(const Version&); ///< Not implemented

}; // Version

#endif // pylith_utils_version_hh


// End of file 
