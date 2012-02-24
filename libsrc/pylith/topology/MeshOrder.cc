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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "MeshOrder.hh" // implementation of class methods

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
ALE::MeshOrder::MeshOrder(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
ALE::MeshOrder::~MeshOrder(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Determine order from pre-existing mesh.
void
ALE::MeshOrder::initialize(const ALE::Obj<mesh_type>& mesh)
{ // initialize
  assert(!mesh.isNull());

  const ALE::Obj<mesh_type::label_sequence>& cells = mesh->heightStratum(0);
  assert(!cells.isNull());

  const Obj<mesh_type::label_sequence>& vertices = mesh->depthStratum(0);
  assert(!vertices.isNull());

  if (mesh->hasLabel("censored depth")) {
    // Count number of cells in censored depth (normal cells).
    const Obj<mesh_type::label_sequence>& cellsNormal = mesh->getLabelStratum("censored depth", mesh->depth());
    assert(!cellsNormal.isNull());
    const mesh_type::label_sequence::iterator cellsNormalEnd = cellsNormal->end();
    const int numCellsNormal = cellsNormal->size();
	
    // Count number of remaining cells (censored cells).
    const int numCellsCensored = cells->size() - numCellsNormal;
	
    // Get number of normal vertices.
    assert(!mesh->getFactory().isNull());
    Obj<mesh_type::numbering_type> vNumbering = mesh->getFactory()->getNumbering(mesh, "censored depth", 0);
    assert(!vNumbering.isNull());
    const int numVerticesNormal = vNumbering->size();
	
    // Count number of remaining vertices (censored vertices).
    const int numVerticesCensored = vertices->size() - numVerticesNormal;

    _cellsNormal = ALE::Interval<point_type>(0, numCellsNormal);
    _verticesNormal = ALE::Interval<point_type>(numCellsNormal, numCellsNormal+numVerticesNormal);
    _verticesCensored = ALE::Interval<point_type>(numCellsNormal+numVerticesNormal, numCellsNormal+numVerticesNormal+numVerticesCensored);
    _cellsCensored = ALE::Interval<point_type>(numCellsNormal+numVerticesNormal+numVerticesCensored,
					       numCellsNormal+numVerticesNormal+numVerticesCensored+numCellsCensored);
  } else {
    const int numCells = cells->size();    
    const int numVertices = vertices->size();

    _cellsNormal = ALE::Interval<point_type>(0, numCells);
    _verticesNormal = ALE::Interval<point_type>(numCells, numCells+numVertices);
    _verticesCensored = ALE::Interval<point_type>(numCells+numVertices, numCells+numVertices);
    _cellsCensored = ALE::Interval<point_type>(numCells+numVertices, numCells+numVertices);
  } // if/else
} // initialize

// ----------------------------------------------------------------------
// Set range for normal cells.
void
ALE::MeshOrder::cellsNormal(const point_type min, const point_type max)
{ // cellsNormal
  _cellsNormal = ALE::Interval<point_type>(min, max);
} // cellsNormal

// ----------------------------------------------------------------------
// Get range for normal cells.
const ALE::Interval<ALE::MeshOrder::point_type>&
ALE::MeshOrder::cellsNormal(void) const
{ // cellsNormal
  return _cellsNormal;
} // cellsNormal

// ----------------------------------------------------------------------
// Set range for normal vertices.
void
ALE::MeshOrder::verticesNormal(const point_type min, const point_type max)
{ // verticesNormal
  _verticesNormal = ALE::Interval<point_type>(min, max);
} // verticesNormal

// ----------------------------------------------------------------------
// Get range for normal vertices.
const ALE::Interval<ALE::MeshOrder::point_type>&
ALE::MeshOrder::verticesNormal(void) const
{ // verticesNormal
  return _verticesNormal;
} // verticesNormal

// ----------------------------------------------------------------------
// Set range for censored cells.
void
ALE::MeshOrder::cellsCensored(const point_type min, const point_type max)
{ // cellsCensored
  _cellsCensored = ALE::Interval<point_type>(min, max);
} // cellsCensored

// ----------------------------------------------------------------------
// Get range for censored cells.
const ALE::Interval<ALE::MeshOrder::point_type>&
ALE::MeshOrder::cellsCensored(void) const
{ // cellsCensored
  return _cellsCensored;
} // cellsCensored

// ----------------------------------------------------------------------
// Set range for censored vertices.
void
ALE::MeshOrder::verticesCensored(const point_type min, const point_type max)
{ // verticesCensored
  _verticesCensored = ALE::Interval<point_type>(min, max);
} // verticesCensored

// ----------------------------------------------------------------------
// Get range for censored vertices.
const ALE::Interval<ALE::MeshOrder::point_type>&
ALE::MeshOrder::verticesCensored(void) const
{ // verticesCensored
  return _verticesCensored;
} // verticesCensored


// End of file 
