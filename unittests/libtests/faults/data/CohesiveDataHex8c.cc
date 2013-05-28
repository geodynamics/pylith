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
 * Cells are 0-1 and vertices are 2-13.
 *
 *       2,3,4,5 -------- 6,7,8,9 -------- 10,11,12,13
 *
 *                        ^^^^^^^ Vertices forming fault
 *
 * After adding cohesive elements
 *
 * Cells are 0-1,2 and vertices are 3-18.
 *
 *       3,4,5,6 -------- 7,8,9,10 -- 15,16,17,18 -------- 11,12,13,14
 *
 *                        ^^^^^^^^^^^^^^^^^^^^^^ Cohesive element
 *
 */

#include "CohesiveDataHex8c.hh"

const int pylith::faults::CohesiveDataHex8c::_numVertices = 16;

const int pylith::faults::CohesiveDataHex8c::_spaceDim = 3;

const int pylith::faults::CohesiveDataHex8c::_numCells = 3;

const int pylith::faults::CohesiveDataHex8c::_cellDim = 3;

const PylithScalar pylith::faults::CohesiveDataHex8c::_vertices[] = {
  -2.0, -1.0, -1.0,
  -2.0,  1.0, -1.0,
  -2.0, -1.0,  1.0,
  -2.0,  1.0,  1.0,
   0.0, -1.0, -1.0,
   0.0,  1.0, -1.0,
   0.0, -1.0,  1.0,
   0.0,  1.0,  1.0,
   2.0, -1.0, -1.0,
   2.0,  1.0, -1.0,
   2.0, -1.0,  1.0,
   2.0,  1.0,  1.0,
   0.0, -1.0, -1.0,
   0.0,  1.0, -1.0,
   0.0, -1.0,  1.0,
   0.0,  1.0,  1.0,
};

const int pylith::faults::CohesiveDataHex8c::_numCorners[] = {
  8,
  8,
  8
};

const int pylith::faults::CohesiveDataHex8c::_cells[] = {
  7,  9, 10,  8,  3,  5,  6,  4,
 11, 13, 14, 12, 15, 17, 18, 16,
  8,  7,  9, 10, 16, 18, 17, 15,
};

const int pylith::faults::CohesiveDataHex8c::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataHex8c::_numGroups = 2;

const int pylith::faults::CohesiveDataHex8c::_groupSizes[] = 
  { 8, 8 };

const int pylith::faults::CohesiveDataHex8c::_groups[] = {
  5, 6, 9, 10, 13, 14, 17, 18,
  7, 8, 9, 10, 15, 16, 17, 18
};

const char* pylith::faults::CohesiveDataHex8c::_groupNames[] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataHex8c::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataHex8c::_filename = 
  "data/hex8c.mesh";

pylith::faults::CohesiveDataHex8c::CohesiveDataHex8c(void)
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

pylith::faults::CohesiveDataHex8c::~CohesiveDataHex8c(void)
{}


// End of file
