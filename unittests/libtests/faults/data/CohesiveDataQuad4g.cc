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

const int pylith::faults::CohesiveDataQuad4g::_numCorners[7] = {
  4,
  4,
  4,
  4,
  4,
  4,
  4,
};

const int pylith::faults::CohesiveDataQuad4g::_materialIds[7] = {
  0,  0,  0,  2,  2,
  1,  1,
};

const int pylith::faults::CohesiveDataQuad4g::_numGroups = 2;

const int pylith::faults::CohesiveDataQuad4g::_groupSizes[2] = {
  3+2, 6+4 // vertices+edges
};

const char* pylith::faults::CohesiveDataQuad4g::_groupNames[2] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataQuad4g::_groupTypes[2] = {
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
  numCorners = const_cast<int*>(_numCorners);
  materialIds = const_cast<int*>(_materialIds);
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
  filename = const_cast<char*>(_filename);
} // constructor

pylith::faults::CohesiveDataQuad4g::~CohesiveDataQuad4g(void)
{}


// End of file
