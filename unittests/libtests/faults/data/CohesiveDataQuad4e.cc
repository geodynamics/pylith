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
 * Cells are 0-3, vertices are 4-12.
 *
 *      10 --------11 --------12
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       5 -------- 7 -------- 9
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       4 -------- 6 -------- 8
 *
 * After adding cohesive elements
 *
 * Cells are 0-3,4-5 vertices are 6-17.
 *
 *      12 --------17--13 --------14
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       7 --------16-- 9 --------11
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       6 --------15-- 8 --------10
 */

#include "CohesiveDataQuad4e.hh"

const int pylith::faults::CohesiveDataQuad4e::_numVertices = 12;

const int pylith::faults::CohesiveDataQuad4e::_spaceDim = 2;

const int pylith::faults::CohesiveDataQuad4e::_numCells = 6;

const int pylith::faults::CohesiveDataQuad4e::_cellDim = 2;

const int pylith::faults::CohesiveDataQuad4e::_numCorners[6] = {
  4,
  4,
  4,
  4,
  4,
  4,
};

const int pylith::faults::CohesiveDataQuad4e::_materialIds[6] = {
  0,  0,  0, 0,
  1,  1,
};

const int pylith::faults::CohesiveDataQuad4e::_numGroups = 2;

const int pylith::faults::CohesiveDataQuad4e::_groupSizes[2] = {
  4+2, 6+4 // vertices+edges
};

const char* pylith::faults::CohesiveDataQuad4e::_groupNames[2] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataQuad4e::_groupTypes[2] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataQuad4e::_filename = 
  "data/quad4e.mesh";

pylith::faults::CohesiveDataQuad4e::CohesiveDataQuad4e(void)
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

pylith::faults::CohesiveDataQuad4e::~CohesiveDataQuad4e(void)
{}


// End of file
