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
 * Cells are 0-1, vertices are 2-6.
 *
 * 2   3,4,5  6
 *
 *     ^^^^^ Face in x-y plane
 *
 * After adding cohesive elements
 *
 * Cells are 0-1,2, vertices are 3-10,11-13.
 *
 * 3   4,5,6  8,9,10   7
 *             11,12,13
 *     ^^^^^^^^^^^^ Cohesive element in x-y plane.
 */

#include "CohesiveDataTet4Lagrange.hh"

const int pylith::faults::CohesiveDataTet4Lagrange::_numVertices = 11;

const int pylith::faults::CohesiveDataTet4Lagrange::_spaceDim = 3;

const int pylith::faults::CohesiveDataTet4Lagrange::_numCells = 3;

const int pylith::faults::CohesiveDataTet4Lagrange::_cellDim = 3;

const PylithScalar pylith::faults::CohesiveDataTet4Lagrange::_vertices[] = {
  -1.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  1.0,
   0.0,  1.0,  0.0,
   1.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  1.0,
   0.0,  1.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  1.0,
   0.0,  1.0,  0.0
};

const int pylith::faults::CohesiveDataTet4Lagrange::_numCorners[] = {
  4,
  4,
  9
};

const int pylith::faults::CohesiveDataTet4Lagrange::_cells[] = {
  4,  5,  6,  3,
  8, 10,  9,  7,
  5,  4,  6,  9,  8,  10,  12, 11, 13
};

const int pylith::faults::CohesiveDataTet4Lagrange::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataTet4Lagrange::_numGroups = 2;

const int pylith::faults::CohesiveDataTet4Lagrange::_groupSizes[] = 
  { 5, 9 };

const int pylith::faults::CohesiveDataTet4Lagrange::_groups[] = {
  3, 5, 6, 9, 10,
  4, 5, 6, 8,  9, 10, 11, 12, 13
};

const char* pylith::faults::CohesiveDataTet4Lagrange::_groupNames[] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataTet4Lagrange::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTet4Lagrange::_filename = 
  "data/tet4.mesh";

pylith::faults::CohesiveDataTet4Lagrange::CohesiveDataTet4Lagrange(void)
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

pylith::faults::CohesiveDataTet4Lagrange::~CohesiveDataTet4Lagrange(void)
{}


// End of file
