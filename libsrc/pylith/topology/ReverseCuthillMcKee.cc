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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "ReverseCuthillMcKee.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh

#include <cassert> // USES assert()
#include <stdexcept> // USES std::exception

// ----------------------------------------------------------------------
// Set vertices and cells in mesh.
void
pylith::topology::ReverseCuthillMcKee::reorder(topology::Mesh* mesh)
{ // reorder
  assert(mesh);

  const int commRank = mesh->commRank();
  if (0 == commRank) {
#if 0
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
#else
    assert(0);
#endif
  } // if    
} // reorder


// End of file 
