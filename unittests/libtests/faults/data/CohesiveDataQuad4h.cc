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
// Copyright (c) 2010-2015 University of California, Davis
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
 * After first fault
 *
 *12 ----13 ----14-28 ----15
 * |      |      |  |      |
 * |  0   |  1   | 9|  2   |
 * |      |      |  |      |
 * |      |      |  |      |
 *16 ----17 ----18-29 ----19
 * |      |      |  |      |
 * |  3   |  4   |10|  5   |
 * |      |      |  |      |
 * |      |      |  |      |
 *20 ----21-----22-30 ----23
 * |      |      |  \--11- |
 * |  6   |  7   |     8   |
 * |      |      |         |
 * |      |      |         |
 *24 ----25 ----26--------27
 *
 * After second fault
 *
 * Cells are 0-8,9-11 vertices are 14-34.
 *
 *14 ----15 ----16-30 ----17
 * |      |      |  |      |
 * |  0   |  1   | 9|  2   |
 * |      |      |  |      |
 * |      |      |  |      |
 *18 ----19 ----20-31 ----21
 * |      |      |  |      |
 * |  3   |  4   |10|  5   |
 * |      |      |  |      |
 * |      |      |  |      |
 *33 ----34-----24-32 ----25
 * |  12  | 13 / |  \-11-- |
 *22 ----23---/  |         |
 * |      |   7  |     8   |
 * |  6   |      |         |
 * |      |      |         |
 * |      |      |         |
 *26 ----27 ----28--------29
 *
 */

#include "CohesiveDataQuad4h.hh"

const int pylith::faults::CohesiveDataQuad4h::_numVertices = 21;

const int pylith::faults::CohesiveDataQuad4h::_spaceDim = 2;

const int pylith::faults::CohesiveDataQuad4h::_numCells = 14;

const int pylith::faults::CohesiveDataQuad4h::_cellDim = 2;

const int pylith::faults::CohesiveDataQuad4h::_numCorners[14] = {
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
  4,
  4,
};

const int pylith::faults::CohesiveDataQuad4h::_materialIds[14] = {
  10, 10, 11, 10, 10, 11, 12, 12, 11,
  2, 2, 1, 1, 1, 
};

const int pylith::faults::CohesiveDataQuad4h::_numGroups = 2;

const int pylith::faults::CohesiveDataQuad4h::_groupSizes[2] = {
  4+2, 6+4 // vertices+edges
};

const char* pylith::faults::CohesiveDataQuad4h::_groupNames[2] = {
  "faultB", "faultA"
};

const char* pylith::faults::CohesiveDataQuad4h::_groupTypes[2] = {
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
  numCorners = const_cast<int*>(_numCorners);
  materialIds = const_cast<int*>(_materialIds);
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
  filename = const_cast<char*>(_filename);
} // constructor

pylith::faults::CohesiveDataQuad4h::~CohesiveDataQuad4h(void)
{}


// End of file
