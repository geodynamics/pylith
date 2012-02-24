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
 *
 * Cells are 0-3 and vertices are 4-22.
 *
 *   The fault runs through the center (x = 0)
 *
     1-----7----13     Cells
    /|    /|    /|
   3-|---9-|--15 |      0      2
  /| |  /| |  /| |
 5 | 0-----6----12
 | |/  | |/  | |/
 | 2-----8----14        1      3
 |/    |/    |/
 4----10----16
 *
 * After adding cohesive elements
 *
     1----23----7----13
    /|    /|   /|    /|
   3-|--25-|--9----15 |
  /| |  /| |  | |  /| |
 5 | 0----22----6----12
 | |/  | |/ | |/  | |/
 | 2----24----8----14
 |/    |/   |/    |/
 4----26---10----16
 *
 */

#include "CohesiveDataHex8g.hh"

const int pylith::faults::CohesiveDataHex8g::_numVertices = 24;

const int pylith::faults::CohesiveDataHex8g::_spaceDim = 3;

const int pylith::faults::CohesiveDataHex8g::_numCells = 6;

const int pylith::faults::CohesiveDataHex8g::_cellDim = 3;

const PylithScalar pylith::faults::CohesiveDataHex8g::_vertices[] = {
  -2.0, -1.0, -2.0,
  -2.0,  1.0, -2.0,
  -2.0, -1.0,  0.0,
  -2.0,  1.0,  0.0,
  -2.0, -1.0,  2.0,
  -2.0,  1.0,  2.0,
   0.0, -1.0, -2.0,
   0.0,  1.0, -2.0,
   0.0, -1.0,  0.0,
   0.0,  1.0,  0.0,
   0.0, -1.0,  2.0,
   0.0,  1.0,  2.0,
   2.0, -1.0, -2.0,
   2.0,  1.0, -2.0,
   2.0, -1.0,  0.0,
   2.0,  1.0,  0.0,
   2.0, -1.0,  2.0,
   2.0,  1.0,  2.0,
   0.0, -1.0, -2.0,
   0.0,  1.0, -2.0,
   0.0, -1.0,  0.0,
   0.0,  1.0,  0.0,
   0.0, -1.0,  2.0,
   0.0,  1.0,  2.0,
};

const int pylith::faults::CohesiveDataHex8g::_numCorners[] = {
  8,
  8,
  8,
  8,
  8,
  8
};

const int pylith::faults::CohesiveDataHex8g::_cells[] = {
   4, 22, 23,  5,  6, 24, 25,  7,
   6, 24, 25,  7,  8, 26, 27,  9,
  10, 16, 17, 11, 12, 18, 19, 13,
  12, 18, 19, 13, 14, 20, 21, 15,
  11, 10, 12, 13, 23, 22, 24, 25,
  13, 12, 14, 15, 25, 24, 26, 27,
};

const int pylith::faults::CohesiveDataHex8g::_materialIds[] = {
  0,  0,  0,  0,
  1,  1
};

const int pylith::faults::CohesiveDataHex8g::_numGroups = 2;

const int pylith::faults::CohesiveDataHex8g::_groupSizes[] = 
  { 12, 8 };

const int pylith::faults::CohesiveDataHex8g::_groups[] = {
  10, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 27,
  8, 9, 14, 15, 20, 21, 26, 27
};

const char* pylith::faults::CohesiveDataHex8g::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataHex8g::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataHex8g::_filename = 
  "data/hex8g.mesh";

pylith::faults::CohesiveDataHex8g::CohesiveDataHex8g(void)
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

pylith::faults::CohesiveDataHex8g::~CohesiveDataHex8g(void)
{}


// End of file
