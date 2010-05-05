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

#include <portinfo>

#include "ReverseCuthillMcKee.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh

#include <cassert> // USES assert()
#include <stdexcept> // USES std::exception

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;

// ----------------------------------------------------------------------
// Set vertices and cells in mesh.
void
pylith::topology::ReverseCuthillMcKee::reorder(topology::Mesh* mesh)
{ // reorder
  assert(0 != mesh);

  //logger.stagePush("MeshReordering");

  int rank = 0;
  MPI_Comm_rank(mesh->comm(), &rank);

  if (0 == rank) {
    const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
    assert(!sieveMesh.isNull());
    ALE::Obj<ALE::Ordering<>::perm_type> perm = 
      new ALE::Ordering<>::perm_type(sieveMesh->comm(), sieveMesh->debug());
    ALE::Obj<ALE::Ordering<>::perm_type> reordering = 
      new ALE::Ordering<>::perm_type(sieveMesh->comm(), sieveMesh->debug());
    
    ALE::Ordering<>::calculateMeshReordering(sieveMesh, perm, reordering);
    
    //perm->view("PERMUTATION");
    //reordering->view("REORDERING");
    //sieveMesh->view("MESH BEFORE RELABEL");
    
    sieveMesh->relabel(*reordering);
    //sieveMesh->view("MESH AFTER RELABEL");
  } // if    

    //logger.stagePop();
} // reorder


// End of file 
