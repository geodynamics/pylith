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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/utils/sievefwd.hh
 *
 * @brief Forward declarations for PETSc Sieve objects.
 */

#if !defined(pylith_utils_sievefwd_hh)
#define pylith_utils_sievefwd_hh


/// Namespace for Sieve package.
namespace ALE {
 
  /// PETSc mesh
  template class Mesh<PetscInt,PetscScalar>;
} // ALE

#endif // pylith_utils_sievefwd_hh


// End of file
