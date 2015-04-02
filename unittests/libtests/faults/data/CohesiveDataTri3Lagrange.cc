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
 * Cells are 0-1, 2, vertices are 3-8,9-10.
 *
 *              7 -9- 4
 *             /|     |\
 *            / |     | \
 *           /  |     |  \
 *          /   |     |   \
 *         3    |     |    6
 *          \   |     |   /
 *           \  |     |  /
 *            \ |     | /
 *             \|     |/
 *              8-10- 5
 */

#include "CohesiveDataTri3Lagrange.hh"

const int pylith::faults::CohesiveDataTri3Lagrange::_numVertices = 6;

const int pylith::faults::CohesiveDataTri3Lagrange::_spaceDim = 2;

const int pylith::faults::CohesiveDataTri3Lagrange::_numCells = 3;

const int pylith::faults::CohesiveDataTri3Lagrange::_cellDim = 2;

const int pylith::faults::CohesiveDataTri3Lagrange::_numCorners[3] = {
  3,
  3,
  4,
};

const int pylith::faults::CohesiveDataTri3Lagrange::_materialIds[3] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataTri3Lagrange::_numGroups = 2;

const int pylith::faults::CohesiveDataTri3Lagrange::_groupSizes[2] = {
  5+4, 4+2, // vertices+edges
};

const char* pylith::faults::CohesiveDataTri3Lagrange::_groupNames[2] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataTri3Lagrange::_groupTypes[2] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTri3Lagrange::_filename = 
  "data/tri3.mesh";

pylith::faults::CohesiveDataTri3Lagrange::CohesiveDataTri3Lagrange(void)
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

pylith::faults::CohesiveDataTri3Lagrange::~CohesiveDataTri3Lagrange(void)
{}


// End of file
