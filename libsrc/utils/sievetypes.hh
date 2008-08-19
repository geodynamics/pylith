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

// This turns on memory debugging, see ALE_mem.hh
#define ALE_MEM_LOGGING

/**
 * @file pylith/utils/sievetypes.hh
 *
 * @brief Aliases for Sieve types.
 */

#if !defined(pylith_utils_sievetypes_hh)
#define pylith_utils_sievetypes_hh

#include <petscmesh.hh> // PETSc Mesh

namespace pylith {

  typedef ALE::IMesh Mesh;
  typedef Mesh::sieve_type sieve_type;
  typedef Mesh::real_section_type real_section_type; 
  typedef Mesh::int_section_type int_section_type;
} // pylith

#endif // pylith_utils_sievetypes_hh


// End of file
