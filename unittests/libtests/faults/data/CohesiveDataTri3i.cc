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
 * After adding cohesive elements (THIS SHOULD GENERATE AN
 * ERROR). Using vertices to mark the fault, this test case leads to
 * an ambiguous location for the fault (vertices 4, 5, and 6 are on
 * the fault, so we don't know which sides of triangle 4-5-6 are on
 * the fault).
 *
 * Cells are 0-3, 4-5, vertices are 6-14.
 *
 *        11
 *        / \
 *       /   \
 *      /     \
 *     /       \
 *   10---------- 7 
 *     \       /12 |\
 *      \     / /| | \
 *       \   / / | |  \
 *        \ / /  | |   \
 *         6 /   | |    9
 *          14   | |   /
 *            \  | |  /
 *             \ | | /
 *              \| |/
 *              13-8
 */

#include "CohesiveDataTri3i.hh"

const int pylith::faults::CohesiveDataTri3i::_numVertices = 9;

const int pylith::faults::CohesiveDataTri3i::_spaceDim = 2;

const int pylith::faults::CohesiveDataTri3i::_numCells = 6;

const int pylith::faults::CohesiveDataTri3i::_cellDim = 2;

const int pylith::faults::CohesiveDataTri3i::_numCorners[] = {
  3,
  3,
  3,
  3,
  4,
  4,
};

const int pylith::faults::CohesiveDataTri3i::_materialIds[] = {
  0,  0,  0,  0,
  1,  1
};

const int pylith::faults::CohesiveDataTri3i::_numGroups = 2;

const int pylith::faults::CohesiveDataTri3i::_groupSizes[] = 
  { 5+4, 6+4 }; // vertices+edges

const char* pylith::faults::CohesiveDataTri3i::_groupNames[] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataTri3i::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTri3i::_filename = "data/tri3i.mesh";

pylith::faults::CohesiveDataTri3i::CohesiveDataTri3i(void)
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

pylith::faults::CohesiveDataTri3i::~CohesiveDataTri3i(void)
{}


// End of file
