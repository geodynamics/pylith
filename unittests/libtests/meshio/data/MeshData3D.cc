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

#include "MeshData3D.hh"

const int pylith::meshio::MeshData3D::_numVertices = 14;

const int pylith::meshio::MeshData3D::_spaceDim = 3;

const int pylith::meshio::MeshData3D::_numCells = 2;

const int pylith::meshio::MeshData3D::_cellDim = 3;

const int pylith::meshio::MeshData3D::_numCorners = 8;

const double pylith::meshio::MeshData3D::_vertices[] = {
  -3.0, -1.0,  0.2,
  -3.0, -1.0,  1.3,
  -1.0, -1.2,  0.1,
  -1.0, -1.2,  1.2,
  -3.0,  5.0,  1.3,
  -3.0,  5.0,  0.1,
  -0.5,  4.8,  0.2,
  -0.5,  4.8,  1.4,
   0.5,  7.0,  1.2,
   1.0,  3.1,  1.3,
   3.0,  4.1,  1.4,
   0.5,  7.0, -0.1,
   1.0,  3.0, -0.2,
   3.0,  4.2,  0.1
};

const int pylith::meshio::MeshData3D::_cells[] = {
  6, 12, 13, 11,  7,  9, 10,  8,
  0,  2,  6,  5,  1,  3,  7,  4
};
const int pylith::meshio::MeshData3D::_materialIds[] = {
  1, 0
};

const bool pylith::meshio::MeshData3D::_useIndexZero = true;

pylith::meshio::MeshData3D::MeshData3D(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numCorners = _numCorners;
  vertices = const_cast<double*>(_vertices);
  cells = const_cast<int*>(_cells);
  materialIds = const_cast<int*>(_materialIds);
  useIndexZero = _useIndexZero;
} // constructor

pylith::meshio::MeshData3D::~MeshData3D(void)
{}


// End of file
