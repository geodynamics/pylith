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
 * @file pylith/topology/Mesh.hh
 *
 * @brief Type definitions for use of PETSc mesh in PyLith.
 */

#if !defined(pylith_topology_mesh_hh)
#define pylith_topology_mesh_hh

namespace pylith {
  namespace topology {
    class Mesh;
  } // topology
} // pylith

class pylith::topology::Mesh
{ // Mesh

  // PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :

  typedef ALE::Obj Obj;
  typedef ALE::Field::Mesh Mesh;
  typedef Mesh::point_type point_type;
  typedef Mesh::real_section_type real_section_type;

}; // Mesh


#endif // pylith_topology_mesh_hh


// End of file
