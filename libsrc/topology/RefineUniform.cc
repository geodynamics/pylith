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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "RefineUniform.hh" // implementation of class methods

#include "Mesh.hh" // USES Mesh

#include "CellRefinerTri3.hh" // USES CellRefinerTri3
#include "MeshRefiner.hh" // USES MeshRefiner

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;

// ----------------------------------------------------------------------
// Constructor
pylith::topology::RefineUniform::RefineUniform(void)
{ // constructor
} // constructor
 
// ----------------------------------------------------------------------
// Destructor
pylith::topology::RefineUniform::~RefineUniform(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Refine mesh.
void
pylith::topology::RefineUniform::refine(Mesh* const newMesh,
					const Mesh& mesh,
					const int levels)
{ // refine
  assert(0 != newMesh);

  typedef SieveMesh::point_type point_type;

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh>& newSieveMesh = newMesh->sieveMesh();
  assert(!newSieveMesh.isNull());
  ALE::Obj<SieveMesh::sieve_type> newSieve =
    new SieveMesh::sieve_type(mesh.comm(), mesh.debug());
  newSieveMesh->setSieve(newSieve);

  ALE::CellRefinerTri3 cellSplitter(*sieveMesh);
  ALE::MeshRefiner refiner;
  refiner.refine(newSieveMesh, sieveMesh, cellSplitter);

} // refine
    

// End of file 
