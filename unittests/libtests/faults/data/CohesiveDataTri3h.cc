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

#include "CohesiveDataTri3h.hh"

const int pylith::faults::CohesiveDataTri3h::_numVertices = 11;

const int pylith::faults::CohesiveDataTri3h::_spaceDim = 2;

const int pylith::faults::CohesiveDataTri3h::_numCells = 13;

const int pylith::faults::CohesiveDataTri3h::_cellDim = 2;

const int pylith::faults::CohesiveDataTri3h::_numCorners[13] = {
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3,  3,
  3, 3,
};

const int pylith::faults::CohesiveDataTri3h::_materialIds[13] = {
  0, 0, 0, 0, 
  2, 2, 3, 3, 3, 2, 3,
  1, 1, 
};

const int pylith::faults::CohesiveDataTri3h::_numGroups = 3;

const int pylith::faults::CohesiveDataTri3h::_groupSizes[3] = 
  { 3+2, 2, 4+4 }; // vertices+edges

const char* pylith::faults::CohesiveDataTri3h::_groupNames[3] = {
  "output", "edge", "fault"
};

const char* pylith::faults::CohesiveDataTri3h::_groupTypes[3] = {
  "vertex", "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTri3h::_filename = "data/tri3h.mesh";
const char* pylith::faults::CohesiveDataTri3h::_fault = "fault";
const char* pylith::faults::CohesiveDataTri3h::_edge = "edge";

pylith::faults::CohesiveDataTri3h::CohesiveDataTri3h(void)
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

pylith::faults::CohesiveDataTri3h::~CohesiveDataTri3h(void)
{}


// End of file
