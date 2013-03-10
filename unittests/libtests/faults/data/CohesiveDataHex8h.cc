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
 | 2-----8----14        3      1
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

#include "CohesiveDataHex8h.hh"

const int pylith::faults::CohesiveDataHex8h::_numVertices = 24;

const int pylith::faults::CohesiveDataHex8h::_spaceDim = 3;

const int pylith::faults::CohesiveDataHex8h::_numCells = 6;

const int pylith::faults::CohesiveDataHex8h::_cellDim = 3;

const PylithScalar pylith::faults::CohesiveDataHex8h::_vertices[] = {
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

const int pylith::faults::CohesiveDataHex8h::_numCorners[] = {
  8,
  8,
  8,
  8,
  8,
  8
};

const int pylith::faults::CohesiveDataHex8h::_cells[] = {
   6, 24, 25,  7,  8, 26, 27,  9,
  21, 20, 22, 23, 15, 14, 16, 17,
  13, 19, 21, 15, 12, 18, 20, 14,
  29, 28, 10, 11, 27, 26,  8,  9,
  15, 13, 12, 14, 27, 25, 24, 26,
  17, 15, 14, 16, 29, 27, 26, 28,
};

const int pylith::faults::CohesiveDataHex8h::_materialIds[] = {
  0,  0,  0,  0,
  1,  1
};

const int pylith::faults::CohesiveDataHex8h::_numGroups = 2;

const int pylith::faults::CohesiveDataHex8h::_groupSizes[] = 
  { 8, 12 };

const int pylith::faults::CohesiveDataHex8h::_groups[] = {
  10, 11, 16, 17, 22, 23, 28, 29,
  12, 13, 14, 15, 16, 17, 24, 25, 26, 27, 28, 29
};

const char* pylith::faults::CohesiveDataHex8h::_groupNames[] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataHex8h::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataHex8h::_filename = 
  "data/hex8h.mesh";

pylith::faults::CohesiveDataHex8h::CohesiveDataHex8h(void)
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

pylith::faults::CohesiveDataHex8h::~CohesiveDataHex8h(void)
{}


// End of file
