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
 * Cells are 0-3, vertices are 4-9.
 *
 *       5 -------- 7 -------- 9
 *       |        / |        / |
 *       |       /  |       /  |
 *       |  0   /   |  1   /   |
 *       |     /    |     /    |
 *       |    /     |    /     |
 *       |   /      |   /      |
 *       |  /   3   |  /   2   |
 *       | /        | /        |
 *       4 -------- 6 -------- 8
 *
 * After adding cohesive elements
 *
 * Cells are 0-3, 4, vertices are 5-12.
 *
 *       6 -------- 8 --12 --------10
 *       |        / |    |        / |
 *       |       /  |    |       /  |
 *       |  0   /   |    |  1   /   |
 *       |     /    |  4 |     /    |
 *       |    /     |    |    /     |
 *       |   /      |    |   /      |
 *       |  /   3   |    |  /   2   |
 *       | /        |    | /        |
 *       5 -------- 7 --11--------- 9
 */

#include "CohesiveDataTri3f.hh"

const int pylith::faults::CohesiveDataTri3f::_numVertices = 8;

const int pylith::faults::CohesiveDataTri3f::_spaceDim = 2;

const int pylith::faults::CohesiveDataTri3f::_numCells = 5;

const int pylith::faults::CohesiveDataTri3f::_cellDim = 2;

const PylithScalar pylith::faults::CohesiveDataTri3f::_vertices[] = {
 -2.0, -1.0,
 -2.0,  1.0,
  0.0, -1.0,
  0.0,  1.0,
  2.0, -1.0,
  2.0,  1.0,
  0.0, -1.0,
  0.0,  1.0,
};

const int pylith::faults::CohesiveDataTri3f::_numCorners[] = {
  3,
  3,
  3,
  3,
  4,
};

const int pylith::faults::CohesiveDataTri3f::_cells[] = {
  5,  8,  6,
 10, 12, 11,
 11,  9, 10,
  7,  8,  5,
  7,  8, 11, 12,
};

const int pylith::faults::CohesiveDataTri3f::_materialIds[] = {
  0,  0,  2,  2,
  1,  1
};

const int pylith::faults::CohesiveDataTri3f::_numGroups = 2;

const int pylith::faults::CohesiveDataTri3f::_groupSizes[] = 
  { 4, 4 };

const int pylith::faults::CohesiveDataTri3f::_groups[] = {
  6, 8, 10, 12,
  7, 8, 11, 12
};

const char* pylith::faults::CohesiveDataTri3f::_groupNames[] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataTri3f::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTri3f::_filename = "data/tri3f.mesh";

pylith::faults::CohesiveDataTri3f::CohesiveDataTri3f(void)
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

pylith::faults::CohesiveDataTri3f::~CohesiveDataTri3f(void)
{}


// End of file
