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

/* Original mesh
 *
 *   Cells are 0-1, vertices are 2-4.
 *
 *   2 -------- 3 -------- 4
 *
 * After adding cohesive elements
 *   3 -------- 6 -2- 4 -------- 5
 */

#include "CohesiveDataLine2Lagrange.hh"

const int pylith::faults::CohesiveDataLine2Lagrange::_numVertices = 5;

const int pylith::faults::CohesiveDataLine2Lagrange::_spaceDim = 1;

const int pylith::faults::CohesiveDataLine2Lagrange::_numCells = 3;

const int pylith::faults::CohesiveDataLine2Lagrange::_cellDim = 1;

const PylithScalar pylith::faults::CohesiveDataLine2Lagrange::_vertices[] = {
  -1.0,
   0.0,
   1.0,
   0.0,
   0.0
};

const int pylith::faults::CohesiveDataLine2Lagrange::_numCorners[] = {
  2,
  2,
  3
};

const int pylith::faults::CohesiveDataLine2Lagrange::_cells[] = {
       3,  6,
       4,  5,
       4,  6,  7
};

const int pylith::faults::CohesiveDataLine2Lagrange::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataLine2Lagrange::_numGroups = 2;

const int pylith::faults::CohesiveDataLine2Lagrange::_groupSizes[] = 
  { 3, 3 };

const int pylith::faults::CohesiveDataLine2Lagrange::_groups[] = {
  3, 4, 6,
  4, 6, 7
};

const char* pylith::faults::CohesiveDataLine2Lagrange::_groupNames[] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataLine2Lagrange::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataLine2Lagrange::_filename = 
  "data/line2.mesh";

pylith::faults::CohesiveDataLine2Lagrange::CohesiveDataLine2Lagrange(void)
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

pylith::faults::CohesiveDataLine2Lagrange::~CohesiveDataLine2Lagrange(void)
{}


// End of file
