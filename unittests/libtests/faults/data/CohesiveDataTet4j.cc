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
 * Cells are 0-5,6 vertices are 7-18.
 */

#include "CohesiveDataTet4j.hh"

const int pylith::faults::CohesiveDataTet4j::_numVertices = 12;

const int pylith::faults::CohesiveDataTet4j::_spaceDim = 3;

const int pylith::faults::CohesiveDataTet4j::_numCells = 7;

const int pylith::faults::CohesiveDataTet4j::_cellDim = 3;

const PylithScalar pylith::faults::CohesiveDataTet4j::_vertices[] = {
  -2.0, -1.0,  0.0,
  -2.0,  0.0,  0.0,
  -2.0,  0.0,  1.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  0.0,
   0.0,  0.0,  1.0,
   2.0, -1.0,  0.0,
   2.0,  0.0,  0.0,
   2.0,  0.0,  1.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  0.0,
   0.0,  0.0,  1.0,
};

const int pylith::faults::CohesiveDataTet4j::_numCorners[] = {
  4,
  4,
  4,
  4,
  4,
  4,
  6
};

const int pylith::faults::CohesiveDataTet4j::_cells[] = {
   8, 10,  9,  7,
  18, 16, 14, 15,
  17, 14, 18, 16,
  11, 10, 12,  8,
  16, 14, 15, 13,
   8, 12,  9, 10,
  10, 11, 12, 16, 17, 18
};

const int pylith::faults::CohesiveDataTet4j::_materialIds[] = {
  0, 2, 2, 0, 2, 0,
  1, 1
};

const int pylith::faults::CohesiveDataTet4j::_numGroups = 2;

const int pylith::faults::CohesiveDataTet4j::_groupSizes[] = 
  { 8, 6 };

const int pylith::faults::CohesiveDataTet4j::_groups[] = {
  7,  9, 10, 12, 13, 15, 16, 18,
 10, 11, 12, 16, 17, 18
};

const char* pylith::faults::CohesiveDataTet4j::_groupNames[] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataTet4j::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTet4j::_filename = "data/tet4j.mesh";

pylith::faults::CohesiveDataTet4j::CohesiveDataTet4j(void)
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

pylith::faults::CohesiveDataTet4j::~CohesiveDataTet4j(void)
{}


// End of file
