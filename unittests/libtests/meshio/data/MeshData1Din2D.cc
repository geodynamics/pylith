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

#include "MeshData1Din2D.hh"

const int pylith::meshio::MeshData1Din2D::_numVertices = 4;

const int pylith::meshio::MeshData1Din2D::_spaceDim = 2;

const int pylith::meshio::MeshData1Din2D::_numCells = 3;

const int pylith::meshio::MeshData1Din2D::_cellDim = 1;

const int pylith::meshio::MeshData1Din2D::_numCorners = 2;

const double pylith::meshio::MeshData1Din2D::_vertices[] = {
  -3.0, -1.2,
   1.0, -1.0,
   2.6,  3.1,
   1.8, -4.0
};

const int pylith::meshio::MeshData1Din2D::_cells[] = {
       3,  1,
       0,  1,
       1,  2
};
const int pylith::meshio::MeshData1Din2D::_materialIds[] = {
  1, 0, 1
};

const bool pylith::meshio::MeshData1Din2D::_useIndexZero = true;

pylith::meshio::MeshData1Din2D::MeshData1Din2D(void)
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

pylith::meshio::MeshData1Din2D::~MeshData1Din2D(void)
{}


// End of file
