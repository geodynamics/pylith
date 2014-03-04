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
 * Cells are 0-5, vertices are 6-17.
 *
 * 6 ---- 7 ---- 8
 * |      |      |
 * |      |      |
 * |      |      |
 * |      |      |
 * 9 ----10 ----11
 * |      |      |
 * |      |      |
 * |      |      |
 * |      |      |
 *12 ----13 ----14
 * |      |      |
 * |      |      |
 * |      |      |
 * |      |      |
 *15 ----16 ----17
 *
 * After adding fault.
 *
 *
 * Cells are 0-8, vertices are 9-23, edges are 24-47.
 *
 * 9 ----10 ------11
 * |      |       |
 * |  0   |   1   |
 * |      |       |
 * |      |  6    |
 *12 ----13-21----14
 * |      |  |      |
 * |      |7 |      |
 * |  2   |  |  3   |
 * |      |  |      |
 *15 ----16-22----17
 * |      |  |      |
 * |  4   |  |  5   |
 * |      |8 |      |
 * |      |  |      |
 *18 ----19-23 ----20
 *
 */

#include "CohesiveDataQuad4i.hh"

const int pylith::faults::CohesiveDataQuad4i::_numVertices = 15;

const int pylith::faults::CohesiveDataQuad4i::_spaceDim = 2;

const int pylith::faults::CohesiveDataQuad4i::_numCells = 9;

const int pylith::faults::CohesiveDataQuad4i::_cellDim = 2;

const int pylith::faults::CohesiveDataQuad4i::_numCorners[9] = {
  4, 4, 4, 4, 4, 4,
  3, 4, 4,
};

const int pylith::faults::CohesiveDataQuad4i::_materialIds[9] = {
  10, 10, 10, 11, 10, 11,
  1, 1, 1,
};

const int pylith::faults::CohesiveDataQuad4i::_numGroups = 3;

const int pylith::faults::CohesiveDataQuad4i::_groupSizes[3] = {
  3+2, 4+2, 6+4, // vertices+edges
};

const char* pylith::faults::CohesiveDataQuad4i::_groupNames[3] = {
  "output2", "output1", "fault",
};

const char* pylith::faults::CohesiveDataQuad4i::_groupTypes[3] = {
  "vertex", "vertex", "vertex",
};

const char* pylith::faults::CohesiveDataQuad4i::_filename = 
  "data/quad4i.mesh";

pylith::faults::CohesiveDataQuad4i::CohesiveDataQuad4i(void)
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

pylith::faults::CohesiveDataQuad4i::~CohesiveDataQuad4i(void)
{}


// End of file
