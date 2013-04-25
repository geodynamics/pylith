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
 *       4 -------- 6 --10 -------- 8
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       3 -------- 5 -- 9 -------- 7
 */

#include "CohesiveDataQuad4.hh"

const int pylith::faults::CohesiveDataQuad4::_numVertices = 8;

const int pylith::faults::CohesiveDataQuad4::_spaceDim = 2;

const int pylith::faults::CohesiveDataQuad4::_numCells = 3;

const int pylith::faults::CohesiveDataQuad4::_cellDim = 2;

const PylithScalar pylith::faults::CohesiveDataQuad4::_vertices[] = {
  -2.0, -1.0,
  -2.0,  1.0,
   0.0, -1.0,
   0.0,  1.0,
   2.0, -1.0,
   2.0,  1.0,
   0.0, -1.0,
   0.0,  1.0,
};

const int pylith::faults::CohesiveDataQuad4::_numCorners[] = {
  4,
  4,
  4
};

const int pylith::faults::CohesiveDataQuad4::_cells[] = {
  3,  5,  6,  4,
  7,  8, 10,  9,
  5,  6,  9, 10,
};

const int pylith::faults::CohesiveDataQuad4::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataQuad4::_numGroups = 2;

const int pylith::faults::CohesiveDataQuad4::_groupSizes[] = 
  { 4, 4 };

const int pylith::faults::CohesiveDataQuad4::_groups[] = {
  4, 6, 8, 10,
  5, 6, 9, 10
};

const char* pylith::faults::CohesiveDataQuad4::_groupNames[] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataQuad4::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataQuad4::_filename = 
  "data/quad4.mesh";

pylith::faults::CohesiveDataQuad4::CohesiveDataQuad4(void)
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

pylith::faults::CohesiveDataQuad4::~CohesiveDataQuad4(void)
{}


// End of file
