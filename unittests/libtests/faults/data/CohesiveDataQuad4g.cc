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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/* Original mesh
 *
 * Cells are 0-4, vertices are 5-15.
 *
 *       7 --------12 --------15
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |     0    |    2     |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       6 --------11 --------14
 *       |        / |          |
 *       |   3   /  |          |
 *       |      /   |          |
 *       5-----9    |          |
 *             |    |    1     |
 *             | 4  |          |
 *             |    |          |
 *             |    |          |
 *             8---10 --------13
 *
 * After adding cohesive elements
 *
 * Cells are 0-4,19-20 vertices are 5-18.
 *
 *       7 --------12 -- 18 --------15
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       6 --------11 -- 17 --------14
 *       |        / |     |          |
 *       |       /  |     |          |
 *       |      /   |     |          |
 *       5-----9    |     |          |
 *             |    |     |          |
 *             |    |     |          |
 *             |    |     |          |
 *             |    |     |          |
 *             8---10 -- 16 --------13
 */

#include "CohesiveDataQuad4g.hh"

const int pylith::faults::CohesiveDataQuad4g::_numVertices = 14;

const int pylith::faults::CohesiveDataQuad4g::_spaceDim = 2;

const int pylith::faults::CohesiveDataQuad4g::_numCells = 7;

const int pylith::faults::CohesiveDataQuad4g::_cellDim = 2;

const PylithScalar pylith::faults::CohesiveDataQuad4g::_vertices[] = {
  -2.0, -2.0,
  -2.0,  0.0,
  -2.0,  2.0,
  -1.0, -2.0,
  -1.0, -1.0,
   0.0, -2.0,
   0.0,  0.0,
   0.0,  2.0,
   2.0, -2.0,
   2.0,  0.0,
   2.0,  2.0,
   0.0, -2.0,
   0.0,  0.0,
   0.0,  2.0,
};

const int pylith::faults::CohesiveDataQuad4g::_numCorners[] = {
  4,
  4,
  4,
  4,
  4,
  4,
  4,
};

const int pylith::faults::CohesiveDataQuad4g::_cells[] = {
  6, 11, 12,  7,
 17, 16, 13, 14,
 15, 18, 17, 14,
  6,  5,  9, 11,
  9,  8, 10, 11,
 10, 11, 16, 17,
 11, 12, 17, 18,
};

const int pylith::faults::CohesiveDataQuad4g::_materialIds[] = {
  0,  0,  0,  2,  2,
  1,  1,
};

const int pylith::faults::CohesiveDataQuad4g::_numGroups = 2;

const int pylith::faults::CohesiveDataQuad4g::_groupSizes[] = 
  { 6, 3 };

const int pylith::faults::CohesiveDataQuad4g::_groups[] = {
  10, 11, 12, 16, 17, 18,
  13, 14, 15
};

const char* pylith::faults::CohesiveDataQuad4g::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataQuad4g::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataQuad4g::_filename = 
  "data/quad4g.mesh";

pylith::faults::CohesiveDataQuad4g::CohesiveDataQuad4g(void)
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

pylith::faults::CohesiveDataQuad4g::~CohesiveDataQuad4g(void)
{}


// End of file
