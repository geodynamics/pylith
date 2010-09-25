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

  // Set material ids
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = 
    cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = 
    cells->end();
  const ALE::Obj<SieveMesh::label_type>& materials =
    sieveMesh->getLabel("material-id");

  const ALE::Obj<SieveMesh::label_sequence>& newCells = 
    newSieveMesh->heightStratum(0);
  assert(!newCells.isNull());
  const SieveMesh::label_sequence::iterator newCellsBegin = 
    newCells->begin();
  const SieveMesh::label_sequence::iterator newCellsEnd = 
    newCells->end();
  const ALE::Obj<SieveMesh::label_type>& newMaterials =
    newSieveMesh->createLabel("material-id");
  
  for (SieveMesh::label_sequence::const_iterator c_iter = cellsBegin,
	 cNew_iter = newCellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    const int numNewCellsPerCell = cellSplitter.numNewCells(*c_iter);
    const int material = sieveMesh->getValue(materials, *c_iter);
    
    for(int i=0; i < numNewCellsPerCell; ++i, ++cNew_iter)
      newSieveMesh->setValue(newMaterials, *cNew_iter, material);
  } // for
  
} // refine
    

// End of file 
