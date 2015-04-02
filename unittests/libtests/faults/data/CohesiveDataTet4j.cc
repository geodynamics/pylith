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
 * Cells are 0-5,6 vertices are 7-18.
 */

#include "CohesiveDataTet4j.hh"

const int pylith::faults::CohesiveDataTet4j::_numVertices = 12;

const int pylith::faults::CohesiveDataTet4j::_spaceDim = 3;

const int pylith::faults::CohesiveDataTet4j::_numCells = 7;

const int pylith::faults::CohesiveDataTet4j::_cellDim = 3;

const int pylith::faults::CohesiveDataTet4j::_numCorners[7] = {
  4,
  4,
  4,
  4,
  4,
  4,
  6,
};

const int pylith::faults::CohesiveDataTet4j::_materialIds[7] = {
  0, 2, 2, 0, 2, 0,
  1,
};

const int pylith::faults::CohesiveDataTet4j::_numGroups = 2;

const int pylith::faults::CohesiveDataTet4j::_groupSizes[2] = {
  8+10+4, 6+6+2 // vertices+edges+faces
};

const char* pylith::faults::CohesiveDataTet4j::_groupNames[2] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataTet4j::_groupTypes[2] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTet4j::_filename = "data/tet4j.mesh";

pylith::faults::CohesiveDataTet4j::CohesiveDataTet4j(void)
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

pylith::faults::CohesiveDataTet4j::~CohesiveDataTet4j(void)
{}


// End of file
