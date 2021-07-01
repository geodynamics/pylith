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
 * @file libsrc/utils/PetscVersion.hh
 *
 * @brief C++ object for PETSc version information.
 */

#if !defined(pylith_utils_petscversion_hh)
#define pylith_utils_petscversion_hh

// Include directives ---------------------------------------------------
#include "utilsfwd.hh" // forward declarations

// Version ----------------------------------------------------------
/** @brief C++ object for getting version info.
 */
class pylith::utils::PetscVersion
{ // PetscVersion
  friend class TestPetscVersion; // unit testing

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
  
// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  PetscVersion(const PetscVersion&); ///< Not implemented
  const PetscVersion& operator=(const PetscVersion&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  static const bool _isRelease; ///< Is source code from a release?
  static const char* _version; ///< Version number.
  static const char* _gitRevision; ///< GIT revision.
  static const char* _gitDate; ///< Date of GIT revision.
  static const char* _gitBranch; ///< GIT branch.

  static const char* _petscDir; ///< PETSC_DIR
  static const char* _petscArch; ///< PETSC_ARCH
  
}; // PetscVersion

#endif // pylith_utils_petscversion_hh


// End of file 
