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

#include "MeshDataCubitQuad.hh"

const int pylith::meshio::MeshDataCubitQuad::_numVertices = 6;

const int pylith::meshio::MeshDataCubitQuad::_spaceDim = 2;

const int pylith::meshio::MeshDataCubitQuad::_numCells = 2;

const int pylith::meshio::MeshDataCubitQuad::_cellDim = 2;

const int pylith::meshio::MeshDataCubitQuad::_numCorners = 4;

const double pylith::meshio::MeshDataCubitQuad::_vertices[] = {
   0.0,  0.0,
   1.0,  0.0,
   1.0,  1.0,
   0.0,  1.0,
   2.0,  0.0,
   2.0,  1.0
};

const int pylith::meshio::MeshDataCubitQuad::_cells[] = {
  0,  1,  2,  3,
  1,  4,  5,  2
};
const int pylith::meshio::MeshDataCubitQuad::_materialIds[] = {
  10, 11
};

const int pylith::meshio::MeshDataCubitQuad::_numGroups = 2;

const int pylith::meshio::MeshDataCubitQuad::_groupSizes[] = 
  { 2, 2 };

const int pylith::meshio::MeshDataCubitQuad::_groups[] = {
  0, 3,
  4, 5
};

const char* pylith::meshio::MeshDataCubitQuad::_groupNames[] = {
  "100", "101"
};

const char* pylith::meshio::MeshDataCubitQuad::_groupTypes[] = {
  "vertex", "vertex"
};

const bool pylith::meshio::MeshDataCubitQuad::_useIndexZero = true;

pylith::meshio::MeshDataCubitQuad::MeshDataCubitQuad(void)
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

pylith::meshio::MeshDataCubitQuad::~MeshDataCubitQuad(void)
{}


// End of file
