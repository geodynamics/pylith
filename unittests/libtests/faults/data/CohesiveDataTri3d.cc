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
 * Cells are 0-3, vertices are 4-9.
 *
 *
 *         9
 *        / \
 *       /   \
 *      /     \
 *     /       \
 *    8---------5
 *     \       /|\
 *      \     / | \
 *       \   /  |  \
 *        \ /   |   \
 *         4    |    7
 *          \   |   /
 *           \  |  /
 *            \ | /
 *             \|/
 *              6
 *
 *
 * After adding cohesive elements
 *
 * Cells are 0-3, 13-14, vertices are 4-12.
 *
 *         9
 *        / \
 *       /   \
 *      /     \
 *     /       \
 *    8---------  5
 *    |          /|
 *   12--------10 |
 *     \       /| |\
 *      \     / | | \
 *       \   /  | |  \
 *        \ /   | |   \
 *         4    | |    7
 *          \   | |   /
 *           \  | |  /
 *            \ | | /
 *             \| |/
 *             11-6
 */

#include "CohesiveDataTri3d.hh"

const int pylith::faults::CohesiveDataTri3d::_numVertices = 9;

const int pylith::faults::CohesiveDataTri3d::_spaceDim = 2;

const int pylith::faults::CohesiveDataTri3d::_numCells = 6;

const int pylith::faults::CohesiveDataTri3d::_cellDim = 2;

const double pylith::faults::CohesiveDataTri3d::_vertices[] = {
 -1.0,  0.0,
  0.0,  1.0,
  0.0, -1.0,
  1.0,  0.0,
 -2.0,  1.0,
 -1.0,  2.0,
  0.0,  1.0,
  0.0, -1.0,
 -2.0,  1.0,
};

const int pylith::faults::CohesiveDataTri3d::_numCorners[] = {
  3,
  3,
  3,
  3,
  4,
  4,
};

const int pylith::faults::CohesiveDataTri3d::_cells[] = {
  4, 11, 10,
  5,  6,  7,
 12,  4, 10,
  8,  5,  9,
  8,  5, 12, 10,
  5,  6, 10, 11,
};

const int pylith::faults::CohesiveDataTri3d::_materialIds[] = {
  0,  0,  0,  0,
  1,  1
};

const int pylith::faults::CohesiveDataTri3d::_numGroups = 2;

const int pylith::faults::CohesiveDataTri3d::_groupSizes[] = 
  { 6, 5 };

const int pylith::faults::CohesiveDataTri3d::_groups[] = {
  5, 6, 8, 10, 11, 12,
  5, 6, 7, 10, 11
};

const char* pylith::faults::CohesiveDataTri3d::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataTri3d::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTri3d::_filename = "data/tri3d.mesh";

pylith::faults::CohesiveDataTri3d::CohesiveDataTri3d(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  vertices = const_cast<double*>(_vertices);
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

pylith::faults::CohesiveDataTri3d::~CohesiveDataTri3d(void)
{}


// End of file
