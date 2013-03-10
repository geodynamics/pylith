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
 * Cells are 0-1, vertices are 2-5.
 *
 *              3
 *             /|\
 *            / | \
 *           /  |  \
 *          /   |   \
 *         2    |    5
 *          \   |   /
 *           \  |  /
 *            \ | /
 *             \|/
 *              4
 *
 *
 * After adding cohesive elements
 *
 * Cells are 0-1, 2, vertices are 3-8.
 *
 *              4 -- 7
 *             /|    |\
 *            / |    | \
 *           /  |    |  \
 *          /   |    |   \
 *         3    |    |    6
 *          \   |    |   /
 *           \  |    |  /
 *            \ |    | /
 *             \|    |/
 *              5 -- 8
 */

#include "CohesiveDataTri3b.hh"

const int pylith::faults::CohesiveDataTri3b::_numVertices = 6;

const int pylith::faults::CohesiveDataTri3b::_spaceDim = 2;

const int pylith::faults::CohesiveDataTri3b::_numCells = 3;

const int pylith::faults::CohesiveDataTri3b::_cellDim = 2;

const PylithScalar pylith::faults::CohesiveDataTri3b::_vertices[] = {
 -1.0,  0.0,
  0.0,  1.0,
  0.0, -1.0,
  1.0,  0.0,
  0.0,  1.0,
  0.0, -1.0
};

const int pylith::faults::CohesiveDataTri3b::_numCorners[] = {
  3,
  3,
  4
};

const int pylith::faults::CohesiveDataTri3b::_cells[] = {
  5,  4,  3,
  8,  6,  7,
  5,  4,  8, 7
};

const int pylith::faults::CohesiveDataTri3b::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataTri3b::_numGroups = 2;

const int pylith::faults::CohesiveDataTri3b::_groupSizes[] = 
  { 5, 4 };

const int pylith::faults::CohesiveDataTri3b::_groups[] = {
  4, 5, 6, 7, 8,
  4, 5, 7, 8
};

const char* pylith::faults::CohesiveDataTri3b::_groupNames[] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataTri3b::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTri3b::_filename = "data/tri3b.mesh";

pylith::faults::CohesiveDataTri3b::CohesiveDataTri3b(void)
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

pylith::faults::CohesiveDataTri3b::~CohesiveDataTri3b(void)
{}


// End of file
