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

#include "MeshData2Din3D.hh"

const int pylith::meshio::MeshData2Din3D::_numVertices = 9;

const int pylith::meshio::MeshData2Din3D::_spaceDim = 2;

const int pylith::meshio::MeshData2Din3D::_numCells = 3;

const int pylith::meshio::MeshData2Din3D::_cellDim = 2;

const int pylith::meshio::MeshData2Din3D::_numCorners = 4;

const double pylith::meshio::MeshData2Din3D::_vertices[] = {
  -1.0,  3.0,  0.2,
   1.0,  3.3,  0.5,
  -1.2,  0.9,  0.3,
   0.9,  1.0,  0.4,
   3.0,  2.9, -0.1,
   6.0,  1.2, -0.2,
   3.4, -0.2,  0.1,
   0.1, -1.1,  0.9,
   2.9, -3.1,  0.8
};

const int pylith::meshio::MeshData2Din3D::_cells[] = {
  0,  2,  3,  1,
  4,  3,  6,  5,
  3,  7,  8,  6
};
const int pylith::meshio::MeshData2Din3D::_materialIds[] = {
  0, 1, 0
};

const bool pylith::meshio::MeshData2Din3D::_useIndexZero = true;

pylith::meshio::MeshData2Din3D::MeshData2Din3D(void)
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

pylith::meshio::MeshData2Din3D::~MeshData2Din3D(void)
{}


// End of file
