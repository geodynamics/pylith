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

const int pylith::faults::CohesiveDataTri3g::_numVertices = 9;

const int pylith::faults::CohesiveDataTri3g::_spaceDim = 2;

const int pylith::faults::CohesiveDataTri3g::_numCells = 8;

const int pylith::faults::CohesiveDataTri3g::_cellDim = 2;

const int pylith::faults::CohesiveDataTri3g::_numCorners[8] = {
  3, 3, 3, 3, 3, 3,
  3, 4,
};

const int pylith::faults::CohesiveDataTri3g::_materialIds[8] = {
  0,  2,  0,  2,  0,  2,
  1, 1,
};

const int pylith::faults::CohesiveDataTri3g::_numGroups = 2;

const int pylith::faults::CohesiveDataTri3g::_groupSizes[2] = 
  { 5+3, 4+2 }; // vertices+edges

const char* pylith::faults::CohesiveDataTri3g::_groupNames[2] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataTri3g::_groupTypes[2] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTri3g::_filename = "data/tri3g.mesh";

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
} // constructor

pylith::faults::CohesiveDataTri3g::~CohesiveDataTri3g(void)
{}


// End of file
