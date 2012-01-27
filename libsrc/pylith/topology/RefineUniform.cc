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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "RefineUniform.hh" // implementation of class methods

#include "Mesh.hh" // USES Mesh

#include "CellRefinerTri3.hh" // USES CellRefinerTri3
#include "CellRefinerQuad4.hh" // USES CellRefinerQuad4
#include "CellRefinerTet4.hh" // USES CellRefinerTet4
#include "CellRefinerHex8.hh" // USES CellRefinerHex8
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

  const Obj<SieveMesh::label_sequence>& cells = sieveMesh->heightStratum(0);
  assert(!cells.isNull());
  assert(cells->size() > 0);

  const int numCorners = sieveMesh->getNumCellCorners();
  const int dim = mesh.dimension();

  switch (dim) {
  case 0:
  case 1:
    throw std::runtime_error("Uniform refinement not implemented.");
    break;

  case 2:
    switch (numCorners) {
    case 3: {
      ALE::CellRefinerTri3 cellSplitter(*sieveMesh);
      ALE::MeshRefiner<ALE::CellRefinerTri3> refinement;
      refinement.refine(newSieveMesh, sieveMesh, cellSplitter);
      break;
    } // case 3
    case 4: {
      ALE::CellRefinerQuad4 cellSplitter(*sieveMesh);
      ALE::MeshRefiner<ALE::CellRefinerQuad4> refinement;
      refinement.refine(newSieveMesh, sieveMesh, cellSplitter);
      break;
    } // case 4
    default :
      throw std::runtime_error("Unknown number of corners for cells.");
    } // switch
    break;
    
  case 3:
    switch (numCorners) {
    case 4: {
      ALE::CellRefinerTet4 cellSplitter(*sieveMesh);
      ALE::MeshRefiner<ALE::CellRefinerTet4> refinement;
      refinement.refine(newSieveMesh, sieveMesh, cellSplitter);
      break;
    } // case 4
    case 8: {
      ALE::CellRefinerHex8 cellSplitter(*sieveMesh);
      ALE::MeshRefiner<ALE::CellRefinerHex8> refinement;
      refinement.refine(newSieveMesh, sieveMesh, cellSplitter);
      break;
    } // case 8
    default :
      assert(0);
      throw std::runtime_error("Unknown number of corners for cells.");
    } // switch
    break;

  default :
    assert(0);
    throw std::logic_error("Unknown dimension.");
  } // switch

  // newMesh->view("REFINED MESH");
} // refine
    

// End of file 
