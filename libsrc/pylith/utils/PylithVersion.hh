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
 * @file libsrc/utils/PylithVersion.hh
 *
 * @brief C++ object for PyLith version information.
 */

#if !defined(pylith_utils_pylithversion_hh)
#define pylith_utils_pylithversion_hh

// Include directives ---------------------------------------------------
#include "utilsfwd.hh" // forward declarations

// Version ----------------------------------------------------------
/** @brief C++ object for getting version info.
 */
class pylith::utils::PylithVersion
{ // PylithVersion
  friend class TestPylithVersion; // unit testing

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
  bool isRelease(void);

  /** Get version number.
   *
   * @returns Version number.
   */
  static
  const char* version(void);

  /** Get DOI.
   *
   * @returns DOI.
   */
  static
  const char* doi(void);

  /** Get GIT revision.
   *
   * @returns GIT revision.
   */
  static
  const char* gitRevision(void);

  /** Get GIT hash.
   *
   * @returns GIT hash.
   */
  static
  const char* gitHash(void);

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
  
// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  PylithVersion(const PylithVersion&); ///< Not implemented
  const PylithVersion& operator=(const PylithVersion&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  static const bool _isRelease; ///< Is source code from a release?
  static const char* _version; ///< Version number.
  static const char* _doi; ///< DOI..
  static const char* _gitRevision; ///< GIT revision.
  static const char* _gitDate; ///< Date of GIT revision.
  static const char* _gitHash; ///< GIT hash.
  static const char* _gitBranch; ///< GIT branch.

}; // PylithVersion

#endif // pylith_utils_pylithversion_hh


// End of file 
