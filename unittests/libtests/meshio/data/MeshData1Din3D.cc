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

#include "MeshData1Din3D.hh"

const int pylith::meshio::MeshData1Din3D::_numVertices = 4;

const int pylith::meshio::MeshData1Din3D::_spaceDim = 3;

const int pylith::meshio::MeshData1Din3D::_numCells = 3;

const int pylith::meshio::MeshData1Din3D::_cellDim = 1;

const int pylith::meshio::MeshData1Din3D::_numCorners = 2;

const double pylith::meshio::MeshData1Din3D::_vertices[] = {
  -3.0, -1.2,  0.3,
   1.0, -1.0,  0.0,
   2.6,  3.1, -0.5,
   1.8, -4.0,  1.0
};

const int pylith::meshio::MeshData1Din3D::_cells[] = {
       3,  1,
       0,  1,
       1,  2
};
const int pylith::meshio::MeshData1Din3D::_materialIds[] = {
  1, 1, 0
};

const int pylith::meshio::MeshData1Din3D::_numGroups = 2;

const int pylith::meshio::MeshData1Din3D::_groupSizes[] = 
  { 1, 1 };

const int pylith::meshio::MeshData1Din3D::_groups[] = {
  2,
  1
};

const char* pylith::meshio::MeshData1Din3D::_groupNames[] = {
  "group A", "group B"
};

const char* pylith::meshio::MeshData1Din3D::_groupTypes[] = {
  "vertex", "cell"
};

const bool pylith::meshio::MeshData1Din3D::_useIndexZero = true;

pylith::meshio::MeshData1Din3D::MeshData1Din3D(void)
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

pylith::meshio::MeshData1Din3D::~MeshData1Din3D(void)
{}


// End of file
