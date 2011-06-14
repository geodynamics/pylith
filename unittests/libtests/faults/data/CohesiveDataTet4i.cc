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
 * 4   5,6,7,8  9
 *
 *     ^^^^^ Face in x-y plane
 *
 * After adding cohesive elements
 *
 * Cells are 0-3,14-15, vertices are 4-13.
 *
 * 4   5,6,7,8  10,11,12,13    9
 *
 *     ^^^^^^^^^^^^ Cohesive element in x-y plane.
 */

#include "CohesiveDataTet4i.hh"

const int pylith::faults::CohesiveDataTet4i::_numVertices = 10;

const int pylith::faults::CohesiveDataTet4i::_spaceDim = 3;

const int pylith::faults::CohesiveDataTet4i::_numCells = 6;

const int pylith::faults::CohesiveDataTet4i::_cellDim = 3;

const double pylith::faults::CohesiveDataTet4i::_vertices[] = {
  -1.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  1.0,
   0.0,  1.0,  0.0,
   0.0,  0.0, -1.0,
   1.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  1.0,
   0.0,  1.0,  0.0,
   0.0,  0.0, -1.0,
};

const int pylith::faults::CohesiveDataTet4i::_numCorners[] = {
  4,
  4,
  4,
  4,
  6,
  6,
};

const int pylith::faults::CohesiveDataTet4i::_cells[] = {
  5,  6,  7,  4,
 10, 12, 11,  9,
 13,  9, 12, 10,
  4,  5,  8,  7,
  8,  7,  5, 13, 12, 10,
  7,  6,  5, 12, 11, 10,
};

const int pylith::faults::CohesiveDataTet4i::_materialIds[] = {
  0, 0, 0, 0,
  1, 1
};

const int pylith::faults::CohesiveDataTet4i::_numGroups = 2;

const int pylith::faults::CohesiveDataTet4i::_groupSizes[] = 
  { 8, 2 };

const int pylith::faults::CohesiveDataTet4i::_groups[] = {
  5, 6, 7, 8, 10, 11, 12, 13,
  4, 9
};

const char* pylith::faults::CohesiveDataTet4i::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataTet4i::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTet4i::_filename = "data/tet4i.mesh";

pylith::faults::CohesiveDataTet4i::CohesiveDataTet4i(void)
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

pylith::faults::CohesiveDataTet4i::~CohesiveDataTet4i(void)
{}


// End of file
