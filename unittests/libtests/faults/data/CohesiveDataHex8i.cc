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

/* Original mesh
 *
 * Cells are 0-6 and vertices are 8-43
 */

#include "CohesiveDataHex8i.hh"

const int pylith::faults::CohesiveDataHex8i::_numVertices = 36;

const int pylith::faults::CohesiveDataHex8i::_spaceDim = 3;

const int pylith::faults::CohesiveDataHex8i::_numCells = 9;

const int pylith::faults::CohesiveDataHex8i::_cellDim = 3;

const PylithScalar pylith::faults::CohesiveDataHex8i::_vertices[] = {
  -2.0, -2.0, -2.0,
  -2.0, -1.0, -2.0,
  -3.0,  0.0, -2.0,
  -2.0,  1.0, -2.0,
  -2.0,  2.0, -2.0,
  -2.0, -2.0,  0.0,
  -2.0, -1.0,  0.0,
  -3.0,  0.0,  0.0,
  -2.0,  1.0,  0.0,
  -2.0,  2.0,  0.0,
  -2.0, -1.0,  2.0,
  -3.0,  0.0,  2.0,
  -2.0,  1.0,  2.0,
   0.0, -2.0, -2.0,
   0.0,  0.0, -2.0,
   0.0,  2.0, -2.0,
   0.0, -2.0,  0.0,
   0.0,  0.0,  0.0,
   0.0,  2.0,  0.0,
   0.0,  0.0,  2.0,
   2.0, -2.0, -2.0,
   2.0, -1.0, -2.0,
   3.0,  0.0, -2.0,
   2.0,  1.0, -2.0,
   2.0,  2.0, -2.0,
   2.0, -2.0,  0.0,
   2.0, -1.0,  0.0,
   3.0,  0.0,  0.0,
   2.0,  1.0,  0.0,
   2.0,  2.0,  0.0,
   0.0, -2.0, -2.0,
   0.0,  0.0, -2.0,
   0.0,  2.0, -2.0,
   0.0, -2.0,  0.0,
   0.0,  0.0,  0.0,
   0.0,  2.0,  0.0,
};

const int pylith::faults::CohesiveDataHex8i::_numCorners[] = {
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
};

const int pylith::faults::CohesiveDataHex8i::_cells[] = {
   10,  9, 22, 23, 15, 14, 25, 26,
   16, 11, 12, 17, 15, 10, 23, 26,
   17, 12, 13, 18, 26, 23, 24, 27,
   32, 40, 30, 31, 37, 43, 35, 36,
   43, 37, 32, 40, 44, 38, 33, 41,
   29, 30, 40, 39, 34, 35, 43, 42,
   17, 16, 15, 26, 21, 20, 19, 28,
   22, 25, 26, 23, 39, 40, 43, 42,
   26, 27, 24, 23, 43, 40, 41, 44,
};

const int pylith::faults::CohesiveDataHex8i::_materialIds[] = {
  0,  0,  0,  0,  0,  0,  2,
  1,  1
};

const int pylith::faults::CohesiveDataHex8i::_numGroups = 2;

const int pylith::faults::CohesiveDataHex8i::_groupSizes[] = 
  { 10, 12 };

const int pylith::faults::CohesiveDataHex8i::_groups[] = {
  29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
  22, 23, 24, 25, 26, 27, 39, 40, 41, 42, 43, 44
};

const char* pylith::faults::CohesiveDataHex8i::_groupNames[] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataHex8i::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataHex8i::_filename = 
  "data/hex8i.mesh";

pylith::faults::CohesiveDataHex8i::CohesiveDataHex8i(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  vertices = const_cast<PylithScalar*>(_vertices);
  numCorners = const_cast<int*>(_numCorners);
  cells = const_cast<int*>(_cells);
  materialIds = const_cast<int*>(_materialIds);
  groups = const_cast<int*>(_groups);
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
  filename = const_cast<char*>(_filename);
} // constructor

pylith::faults::CohesiveDataHex8i::~CohesiveDataHex8i(void)
{}


// End of file
