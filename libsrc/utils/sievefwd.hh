// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/**
 * @file pylith/utils/sievefwd.hh
 *
 * @brief Forward declarations for PETSc Sieve objects.
 */

#if !defined(pylith_utils_sievefwd_hh)
#define pylith_utils_sievefwd_hh

/// Namespace for Sieve package.
namespace ALE {
 
  /// ALE::Obj
  template<class T, typename A> class Obj;

  /// PETSc mesh
  class Mesh;
} // ALE

#endif // pylith_utils_sievefwd_hh


// End of file
