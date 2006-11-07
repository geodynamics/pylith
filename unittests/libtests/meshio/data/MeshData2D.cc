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

#include "MeshData2D.hh"

const int pylith::meshio::MeshData2D::_numVertices = 9;

const int pylith::meshio::MeshData2D::_spaceDim = 2;

const int pylith::meshio::MeshData2D::_numCells = 3;

const int pylith::meshio::MeshData2D::_cellDim = 2;

const int pylith::meshio::MeshData2D::_numCorners = 4;

const double pylith::meshio::MeshData2D::_vertices[] = {
  -1.0,  3.0,
   1.0,  3.3,
  -1.2,  0.9,
   0.9,  1.0,
   3.0,  2.9,
   6.0,  1.2,
   3.4, -0.2,
   0.1, -1.1,
   2.9, -3.1
};

const int pylith::meshio::MeshData2D::_cells[] = {
  0,  2,  3,  1,
  4,  3,  6,  5,
  3,  7,  8,  6
};

const bool pylith::meshio::MeshData2D::_useIndexZero = true;

pylith::meshio::MeshData2D::MeshData2D(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numCorners = _numCorners;
  vertices = const_cast<double*>(_vertices);
  cells = const_cast<int*>(_cells);
  useIndexZero = _useIndexZero;
} // constructor

pylith::meshio::MeshData2D::~MeshData2D(void)
{}


// End of file
