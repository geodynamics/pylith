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
 * @file pylith/utils/sievetypes.hh
 *
 * @brief Aliases for Sieve types.
 */

#if !defined(pylith_utils_sievetypes_hh)
#define pylith_utils_sievetypes_hh

#include <petscmesh.hh> // PETSc Mesh

namespace pylith {

  typedef ALE::IMesh<> Mesh;
  typedef ALE::IMesh<ALE::LabelSifter<int, Mesh::point_type> > SubMesh;

} // pylith

#endif // pylith_utils_sievetypes_hh


// End of file
