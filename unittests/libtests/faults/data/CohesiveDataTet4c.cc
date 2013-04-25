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
 * Cells are 0-1,2, vertices are 3-10.
 *
 * 3   4,5,6  8,9,10   7
 *
 *     ^^^^^^^^^^^^ Cohesive element in x-y plane.
 */

#include "CohesiveDataTet4c.hh"

const int pylith::faults::CohesiveDataTet4c::_numVertices = 8;

const int pylith::faults::CohesiveDataTet4c::_spaceDim = 3;

const int pylith::faults::CohesiveDataTet4c::_numCells = 3;

const int pylith::faults::CohesiveDataTet4c::_cellDim = 3;

const PylithScalar pylith::faults::CohesiveDataTet4c::_vertices[] = {
  -1.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  1.0,
   0.0,  1.0,  0.0,
   1.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  1.0,
   0.0,  1.0,  0.0
};

const int pylith::faults::CohesiveDataTet4c::_numCorners[] = {
  4,
  4,
  6
};

const int pylith::faults::CohesiveDataTet4c::_cells[] = {
  6,  5,  3,  4,
  7,  9, 10,  8,
  5,  4,  6,  9,  8, 10,
};

const int pylith::faults::CohesiveDataTet4c::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataTet4c::_numGroups = 2;

const int pylith::faults::CohesiveDataTet4c::_groupSizes[] = 
  { 5, 6 };

const int pylith::faults::CohesiveDataTet4c::_groups[] = {
  3, 5, 6, 9, 10,
  4, 5, 6, 8,  9, 10
};

const char* pylith::faults::CohesiveDataTet4c::_groupNames[] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataTet4c::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTet4c::_filename = "data/tet4c.mesh";

pylith::faults::CohesiveDataTet4c::CohesiveDataTet4c(void)
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

pylith::faults::CohesiveDataTet4c::~CohesiveDataTet4c(void)
{}


// End of file
