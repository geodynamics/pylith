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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/* Original mesh
 *
 */

#include "CohesiveDataTri3g.hh"

const int pylith::faults::CohesiveDataTri3g::_numVertices = 8;

const int pylith::faults::CohesiveDataTri3g::_spaceDim = 2;

const int pylith::faults::CohesiveDataTri3g::_numCells = 7;

const int pylith::faults::CohesiveDataTri3g::_cellDim = 2;

const int pylith::faults::CohesiveDataTri3g::_numCorners[7] = {
  3, 3, 3, 3, 3, 3,
  3,
};

const int pylith::faults::CohesiveDataTri3g::_materialIds[7] = {
  0,  2,  0,  2,  0,  2,
  1,
};

const int pylith::faults::CohesiveDataTri3g::_numGroups = 3;

const int pylith::faults::CohesiveDataTri3g::_groupSizes[3] = 
  { 5+3, 1, 3+2 }; // vertices+edges

const char* pylith::faults::CohesiveDataTri3g::_groupNames[3] = {
  "output", "edge", "fault"
};

const char* pylith::faults::CohesiveDataTri3g::_groupTypes[3] = {
  "vertex", "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTri3g::_filename = "data/tri3g.mesh";
const char* pylith::faults::CohesiveDataTri3g::_fault = "fault";
const char* pylith::faults::CohesiveDataTri3g::_edge = "edge";

pylith::faults::CohesiveDataTri3g::CohesiveDataTri3g(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numCorners = const_cast<int*>(_numCorners);
  materialIds = const_cast<int*>(_materialIds);
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
  filename = const_cast<char*>(_filename);
  fault = const_cast<char*>(_fault);
  edge = const_cast<char*>(_edge);
} // constructor

pylith::faults::CohesiveDataTri3g::~CohesiveDataTri3g(void)
{}


// End of file
