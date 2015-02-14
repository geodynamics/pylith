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
 * 4   5,6,7,8  9
 *
 *     ^^^^^ Face in x-y plane
 *
 * After adding cohesive elements
 *
 * Cells are 0-3,4-5, vertices are 6-15.
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

const int pylith::faults::CohesiveDataTet4i::_numCorners[6] = {
  4,
  4,
  4,
  4,
  6,
  6,
};

const int pylith::faults::CohesiveDataTet4i::_materialIds[6] = {
  0, 0, 0, 0,
  1, 1
};

const int pylith::faults::CohesiveDataTet4i::_numGroups = 2;

const int pylith::faults::CohesiveDataTet4i::_groupSizes[2] = {
  2+0+0, 8+10+4 // vertices+edges+faces
};

const char* pylith::faults::CohesiveDataTet4i::_groupNames[2] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataTet4i::_groupTypes[2] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTet4i::_filename = "data/tet4i.mesh";

pylith::faults::CohesiveDataTet4i::CohesiveDataTet4i(void)
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

pylith::faults::CohesiveDataTet4i::~CohesiveDataTet4i(void)
{}


// End of file
