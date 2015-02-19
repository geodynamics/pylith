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

#include "CohesiveDataTri3c.hh"

const int pylith::faults::CohesiveDataTri3c::_numVertices = 6;

const int pylith::faults::CohesiveDataTri3c::_spaceDim = 2;

const int pylith::faults::CohesiveDataTri3c::_numCells = 3;

const int pylith::faults::CohesiveDataTri3c::_cellDim = 2;

const int pylith::faults::CohesiveDataTri3c::_numCorners[] = {
  3,
  3,
  4
};

const int pylith::faults::CohesiveDataTri3c::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataTri3c::_numGroups = 2;

const int pylith::faults::CohesiveDataTri3c::_groupSizes[] = 
  { 5+4, 4+2 }; // vertices+edges

const char* pylith::faults::CohesiveDataTri3c::_groupNames[] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataTri3c::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTri3c::_filename = "data/tri3c.mesh";

pylith::faults::CohesiveDataTri3c::CohesiveDataTri3c(void)
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

pylith::faults::CohesiveDataTri3c::~CohesiveDataTri3c(void)
{}


// End of file
