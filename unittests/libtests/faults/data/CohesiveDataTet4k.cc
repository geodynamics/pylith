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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "CohesiveDataTet4k.hh"

const int pylith::faults::CohesiveDataTet4k::_numVertices = 12;

const int pylith::faults::CohesiveDataTet4k::_spaceDim = 3;

const int pylith::faults::CohesiveDataTet4k::_numCells = 7;

const int pylith::faults::CohesiveDataTet4k::_cellDim = 3;

const int pylith::faults::CohesiveDataTet4k::_numCorners[7] = {
  4,
  4,
  4,
  4,
  4,
  4,
  6,
};

const int pylith::faults::CohesiveDataTet4k::_materialIds[7] = {
  0, 2, 2, 0, 2, 0,
  1,
};

const int pylith::faults::CohesiveDataTet4k::_numGroups = 2;

const int pylith::faults::CohesiveDataTet4k::_groupSizes[2] = {
  8+10+4, 6+6+2 // vertices+edges+faces
};

const char* pylith::faults::CohesiveDataTet4k::_groupNames[2] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataTet4k::_groupTypes[2] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTet4k::_filename = "data/tet4k.mesh";
const char* pylith::faults::CohesiveDataTet4k::_fault = "fault";
const char* pylith::faults::CohesiveDataTet4k::_edge = "fault_edge";

pylith::faults::CohesiveDataTet4k::CohesiveDataTet4k(void)
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
  fault = const_cast<char*>(_fault);
  edge = const_cast<char*>(_edge);
} // constructor

pylith::faults::CohesiveDataTet4k::~CohesiveDataTet4k(void)
{}


// End of file
