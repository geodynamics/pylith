// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/* Original mesh
 *
 *   0 -------- 1 -------- 2
 *
 * After adding cohesive elements
 *   0 -------- 1 -- 3 -------- 2
 */

#include "CohesiveData1D.hh"

const int pylith::faults::CohesiveData1D::_numVertices = 3;

const int pylith::faults::CohesiveData1D::_spaceDim = 1;

const int pylith::faults::CohesiveData1D::_numCells = 3;

const int pylith::faults::CohesiveData1D::_cellDim = 1;

const double pylith::faults::CohesiveData1D::_vertices[] = {
  -1.0,
   0.0,
   1.0
};

const int pylith::faults::CohesiveData1D::_numCorners[] = {
  2,
  2,
  2
};

const int pylith::faults::CohesiveData1D::_cells[] = {
       0,  1,
       3,  2,
       1,  3,
};

const int pylith::faults::CohesiveData1D::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveData1D::_numGroups = 2;

const int pylith::faults::CohesiveData1D::_groupSizes[] = 
  { 1, 2 };

const int pylith::faults::CohesiveData1D::_groups[] = {
  1,
  0, 1
};

const char* pylith::faults::CohesiveData1D::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveData1D::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveData1D::_filename = "data/meshLine.txt";

pylith::faults::CohesiveData1D::CohesiveData1D(void)
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

pylith::faults::CohesiveData1D::~CohesiveData1D(void)
{}


// End of file
