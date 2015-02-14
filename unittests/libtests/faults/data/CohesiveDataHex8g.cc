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
 *
 * Cells are 0-3 and vertices are 4-22.
 *
 *   The fault runs through the center (x = 0)
 *
     1-----7----13     Cells            y
    /|    /|    /|                      |
   3-|---9-|--15 |      0      2        |
  /| |  /| |  /| |                      |________x
 5 | 0-----6----12                     /
 | |/  | |/  | |/                     /
 | 2-----8----14        1      3     z
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

const int pylith::faults::CohesiveDataHex8g::_numCorners[6] = {
  8,
  8,
  8,
  8,
  8,
  8
};

const int pylith::faults::CohesiveDataHex8g::_materialIds[6] = {
  0,  0,  0,  0,
  1,  1
};

const int pylith::faults::CohesiveDataHex8g::_numGroups = 2;

const int pylith::faults::CohesiveDataHex8g::_groupSizes[2] = {
  8+8+2, 12+14+4 // vertices+edges+faces 
};

const char* pylith::faults::CohesiveDataHex8g::_groupNames[2] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataHex8g::_groupTypes[2] = {
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
  numCorners = const_cast<int*>(_numCorners);
  materialIds = const_cast<int*>(_materialIds);
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
  filename = const_cast<char*>(_filename);
} // constructor

pylith::faults::CohesiveDataHex8g::~CohesiveDataHex8g(void)
{}


// End of file
