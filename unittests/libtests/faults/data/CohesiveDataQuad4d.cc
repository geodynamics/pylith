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
 * Cells are 0-1, vertices are 2-7.
 *
 *       3 -------- 5 -------- 7
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       2 -------- 4 -------- 6
 *
 * After adding cohesive elements
 *
 * Cells are 0-1,2 vertices are 3-10.
 *
 *       4 --------10 -- 6 -------- 8
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       3 -------- 9 -- 5 -------- 7
 */

#include "CohesiveDataQuad4d.hh"

const int pylith::faults::CohesiveDataQuad4d::_numVertices = 8;

const int pylith::faults::CohesiveDataQuad4d::_spaceDim = 2;

const int pylith::faults::CohesiveDataQuad4d::_numCells = 3;

const int pylith::faults::CohesiveDataQuad4d::_cellDim = 2;

const int pylith::faults::CohesiveDataQuad4d::_numCorners[3] = {
  4,
  4,
  4
};

const int pylith::faults::CohesiveDataQuad4d::_materialIds[3] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataQuad4d::_numGroups = 2;

const int pylith::faults::CohesiveDataQuad4d::_groupSizes[2] = {
  4+2, 4+2 // vertices+edges
};

const char* pylith::faults::CohesiveDataQuad4d::_groupNames[2] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataQuad4d::_groupTypes[2] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataQuad4d::_filename = 
  "data/quad4d.mesh";

pylith::faults::CohesiveDataQuad4d::CohesiveDataQuad4d(void)
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

pylith::faults::CohesiveDataQuad4d::~CohesiveDataQuad4d(void)
{}


// End of file
