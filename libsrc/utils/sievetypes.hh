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

#if NEWPYLITHMESH // For use with pylith::topology::Mesh
  typedef ALE::IMesh<> SieveMesh;
  typedef SieveMesh::real_section_type MeshRealSection;
  typedef SieveMesh::int_section_type MeshIntSection;

  typedef ALE::IMesh<ALE::LabelSifter<int, SieveMesh::point_type> > SieveSubMesh;
  typedef SieveSubMesh::real_section_type SubMeshRealSection;
  typedef SieveSubMesh::int_section_type SubMeshIntSection;

#else
  typedef ALE::IMesh<> Mesh;
  typedef ALE::IMesh<ALE::LabelSifter<int, Mesh::point_type> > SubMesh;
  typedef Mesh::sieve_type sieve_type;
  typedef Mesh::real_section_type real_section_type; 
  typedef Mesh::int_section_type int_section_type;
#endif

} // pylith

#endif // pylith_utils_sievetypes_hh


// End of file
