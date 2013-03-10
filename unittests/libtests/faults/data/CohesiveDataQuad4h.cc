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
 * Cells are 0-8, vertices are 9-24.
 *
 * 9 ----10 ----11 ----12
 * |      |      |      |
 * |      |      |      |
 * |      |      |      |
 * |      |      |      |
 *13 ----14 ----15 ----16
 * |      |      |      |
 * |      |      |      |
 * |      |      |      |
 * |      |      |      |
 *17 ----18 ----19 ----20
 * |      |      |      |
 * |      |      |      |
 * |      |      |      |
 * |      |      |      |
 *21 ----22 ----23 ----24
 *
 * After adding cohesive elements
 *
 * Cells are 0-8,9-11 vertices are 12-32.
 *
 *12 ----13 ----28-14 ----15
 * |      |      |  |      |
 * |  0   |  1   | 9|  2   |
 * |      |      |  |      |
 * |      |      |  |      |
 *16 ----17 ----29-18 ----19
 * |      |      |  |      |
 * |  3   |  4   |10|  5   |
 * |      |      |  |      |
 * |      |      |  |      |
 *31 ----32-----30-22 ----23
 * |  11  |      |         |
 *20 ----21      |         |
 * |      |   7  |     8   |
 * |  6   |      |         |
 * |      |      |         |
 * |      |      |         |
 *24 ----25 ----26--------27
 *
 */

#include "CohesiveDataQuad4h.hh"

const int pylith::faults::CohesiveDataQuad4h::_numVertices = 21;

const int pylith::faults::CohesiveDataQuad4h::_spaceDim = 2;

const int pylith::faults::CohesiveDataQuad4h::_numCells = 12;

const int pylith::faults::CohesiveDataQuad4h::_cellDim = 2;

const PylithScalar pylith::faults::CohesiveDataQuad4h::_vertices[] = {
  -3.0,  3.0,
  -1.0,  3.0,
   1.0,  3.0,
   3.0,  3.0,
  -3.0,  1.0,
  -1.0,  1.0,
   1.0,  1.0,
   3.0,  1.0,
  -3.0, -1.0,
  -1.0, -1.0,
   1.0, -1.0,
   3.0, -1.0,
  -3.0, -3.0,
  -1.0, -3.0,
   1.0, -3.0,
   3.0, -3.0,
   1.0,  3.0,
   1.0,  1.0,
   1.0, -1.0,
  -3.0, -1.0,
  -1.0, -1.0,
};

const int pylith::faults::CohesiveDataQuad4h::_numCorners[] = {
  4,
  4,
  4,
  4,
  4,
  4,
  4,
  4,
  4,
  4,
  4,
  4,
};

const int pylith::faults::CohesiveDataQuad4h::_cells[] = {
  12, 16, 17, 13,
  13, 17, 29, 28,
  14, 18, 19, 15,
  16, 31, 32, 17,
  17, 32, 30, 29,
  18, 22, 23, 19,
  20, 24, 25, 21,
  32, 25, 26, 30,
  30, 26, 27, 23,
  14, 18, 28, 29,
  18, 22, 29, 30,
  21, 20, 32, 31,
};

const int pylith::faults::CohesiveDataQuad4h::_materialIds[] = {
  10, 10, 11, 10, 10, 11, 12, 12, 11,
  1, 1, 2,
};

const int pylith::faults::CohesiveDataQuad4h::_numGroups = 2;

const int pylith::faults::CohesiveDataQuad4h::_groupSizes[] = 
  { 4, 6 };

const int pylith::faults::CohesiveDataQuad4h::_groups[] = {
  20, 21, 31, 32,
  14, 18, 22, 28, 29, 30
};

const char* pylith::faults::CohesiveDataQuad4h::_groupNames[] = {
  "faultB", "faultA"
};

const char* pylith::faults::CohesiveDataQuad4h::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataQuad4h::_filename = 
  "data/quad4h.mesh";

pylith::faults::CohesiveDataQuad4h::CohesiveDataQuad4h(void)
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

pylith::faults::CohesiveDataQuad4h::~CohesiveDataQuad4h(void)
{}


// End of file
