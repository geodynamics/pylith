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

#include "MeshDataCubitTet.hh"

const int pylith::meshio::MeshDataCubitTet::_numVertices = 5;

const int pylith::meshio::MeshDataCubitTet::_spaceDim = 3;

const int pylith::meshio::MeshDataCubitTet::_numCells = 2;

const int pylith::meshio::MeshDataCubitTet::_cellDim = 3;

const int pylith::meshio::MeshDataCubitTet::_numCorners = 4;

const double pylith::meshio::MeshDataCubitTet::_vertices[] = {
  -2.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  1.0,  0.0,
   0.0,  0.0,  2.0,
   2.0,  0.0,  0.0
};

const int pylith::meshio::MeshDataCubitTet::_cells[] = {
  0,  1,  2,  3,
  1,  4,  2,  3
};
const int pylith::meshio::MeshDataCubitTet::_materialIds[] = {
  1, 2
};

const int pylith::meshio::MeshDataCubitTet::_numGroups = 2;

const int pylith::meshio::MeshDataCubitTet::_groupSizes[] = 
  { 3, 4 };

const int pylith::meshio::MeshDataCubitTet::_groups[] = {
  1, 2, 3,
  0, 1, 2, 4
};

const char* pylith::meshio::MeshDataCubitTet::_groupNames[] = {
  "100", "101"
};

const char* pylith::meshio::MeshDataCubitTet::_groupTypes[] = {
  "vertex", "vertex"
};

const bool pylith::meshio::MeshDataCubitTet::_useIndexZero = true;

pylith::meshio::MeshDataCubitTet::MeshDataCubitTet(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numCorners = _numCorners;
  vertices = const_cast<double*>(_vertices);
  cells = const_cast<int*>(_cells);
  materialIds = const_cast<int*>(_materialIds);
  groups = const_cast<int*>(_groups);
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
  useIndexZero = _useIndexZero;
} // constructor

pylith::meshio::MeshDataCubitTet::~MeshDataCubitTet(void)
{}


// End of file
