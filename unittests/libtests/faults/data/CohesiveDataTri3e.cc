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
 *
 *         9
 *        / \
 *       /   \
 *      /  2  \
 *     /       \
 *    8---------5
 *     \       /|\
 *      \  3  / | \
 *       \   /  |  \
 *        \ /   |   \
 *         4  0 | 1  7
 *          \   |   /
 *           \  |  /
 *            \ | /
 *             \|/
 *              6
 *
 *
 * After adding cohesive elements
 *
 * Cells are 0-3, 4-5, vertices are 6-14.
 *
 *        11
 *        / \
 *       /   \
 *      /     \
 *     /       \
 *   10---------  7
 *    |          /|
 *   14--------12 |\
 *     \       /| | \
 *      \     / | |  \
 *       \   /  | |   \
 *        \ /   | |    \
 *         6    | |    9
 *          \   | |   /
 *           \  | |  /
 *            \ | | /
 *             \| |/
 *             13-8
 */

#include "CohesiveDataTri3e.hh"

const int pylith::faults::CohesiveDataTri3e::_numVertices = 9;

const int pylith::faults::CohesiveDataTri3e::_spaceDim = 2;

const int pylith::faults::CohesiveDataTri3e::_numCells = 6;

const int pylith::faults::CohesiveDataTri3e::_cellDim = 2;

const PylithScalar pylith::faults::CohesiveDataTri3e::_vertices[] = {
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

const int pylith::faults::CohesiveDataTri3e::_numCorners[] = {
  3,
  3,
  3,
  3,
  4,
  4,
};

const int pylith::faults::CohesiveDataTri3e::_cells[] = {
  6, 13, 12,
  7,  8,  9,
 10,  7, 11,
 14,  6, 12,
  7,  8, 12, 13,
 10,  7, 14, 12,
};

const int pylith::faults::CohesiveDataTri3e::_materialIds[] = {
  0,  0,  0,  0,
  1,  1
};

const int pylith::faults::CohesiveDataTri3e::_numGroups = 2;

const int pylith::faults::CohesiveDataTri3e::_groupSizes[] = 
  { 5, 6 };

const int pylith::faults::CohesiveDataTri3e::_groups[] = {
  7, 8,  9, 12, 13,
  7, 8, 10, 12, 13, 14
};

const char* pylith::faults::CohesiveDataTri3e::_groupNames[] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataTri3e::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTri3e::_filename = "data/tri3e.mesh";

pylith::faults::CohesiveDataTri3e::CohesiveDataTri3e(void)
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

pylith::faults::CohesiveDataTri3e::~CohesiveDataTri3e(void)
{}


// End of file
