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

#include "MeshData1D.hh"

const int pylith::meshio::MeshData1D::_numVertices = 3;

const int pylith::meshio::MeshData1D::_spaceDim = 1;

const int pylith::meshio::MeshData1D::_numCells = 2;

const int pylith::meshio::MeshData1D::_cellDim = 1;

const int pylith::meshio::MeshData1D::_numCorners = 2;

const double pylith::meshio::MeshData1D::_vertices[] = {
  -1.2,
   2.1,
   0.3
};

const int pylith::meshio::MeshData1D::_cells[] = {
       0,  2,
       2,  1,
};

const int pylith::meshio::MeshData1D::_materialIds[] = {
  2,  1
};

const int pylith::meshio::MeshData1D::_numGroups = 3;

const int pylith::meshio::MeshData1D::_groupSizes[] = 
  { 1, 1, 2 };

const int pylith::meshio::MeshData1D::_groups[] = {
  1,
  0,
  0, 1
};

const char* pylith::meshio::MeshData1D::_groupNames[] = {
  "group A", "group B", "group C"
};

const char* pylith::meshio::MeshData1D::_groupTypes[] = {
  "vertex", "cell", "vertex"
};

const bool pylith::meshio::MeshData1D::_useIndexZero = true;

pylith::meshio::MeshData1D::MeshData1D(void)
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

pylith::meshio::MeshData1D::~MeshData1D(void)
{}


// End of file
