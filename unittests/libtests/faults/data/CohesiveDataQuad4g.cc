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
 * Cells are 0-4,5-6 vertices are 7-20.
 *
 *       9 --------14 -- 20 --------17
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       8 --------13 -- 19 --------16
 *       |        / |     |          |
 *       |       /  |     |          |
 *       |      /   |     |          |
 *       7-----11   |     |          |
 *             |    |     |          |
 *             |    |     |          |
 *             |    |     |          |
 *             |    |     |          |
 *            10---12 -- 18 --------15
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
  8, 13, 14,  9,
 19, 18, 15, 16,
 17, 20, 19, 16,
  8,  7, 11, 13,
 11, 10, 12, 13,
 13, 14, 19, 20,
 12, 13, 18, 19,
};

const int pylith::faults::CohesiveDataQuad4g::_materialIds[] = {
  0,  0,  0,  2,  2,
  1,  1,
};

const int pylith::faults::CohesiveDataQuad4g::_numGroups = 2;

const int pylith::faults::CohesiveDataQuad4g::_groupSizes[] = 
  { 3, 6 };

const int pylith::faults::CohesiveDataQuad4g::_groups[] = {
  15, 16, 17,
  12, 13, 14, 18, 19, 20
};

const char* pylith::faults::CohesiveDataQuad4g::_groupNames[] = {
  "output", "fault"
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
