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

#include "MeshDataLagritTet.hh"

const int pylith::meshio::MeshDataLagritTet::_numVertices = 12;

const int pylith::meshio::MeshDataLagritTet::_spaceDim = 3;

const int pylith::meshio::MeshDataLagritTet::_numCells = 12;

const int pylith::meshio::MeshDataLagritTet::_cellDim = 3;

const int pylith::meshio::MeshDataLagritTet::_numCorners = 4;

const double pylith::meshio::MeshDataLagritTet::_vertices[] = {
  0.00000E+000,  -5.00000E-001,  -5.00000E-001,
  0.00000E+000,  -5.00000E-001,   5.00000E-001,
  1.00000E+000,  -5.00000E-001,  -5.00000E-001,
  1.00000E+000,  -5.00000E-001,   5.00000E-001,
  0.00000E+000,   5.00000E-001,  -5.00000E-001,
  0.00000E+000,   5.00000E-001,   5.00000E-001,
  1.00000E+000,   5.00000E-001,  -5.00000E-001,
  1.00000E+000,   5.00000E-001,   5.00000E-001,
 -1.00000E+000,  -5.00000E-001,  -5.00000E-001,
 -1.00000E+000,  -5.00000E-001,   5.00000E-001,
 -1.00000E+000,   5.00000E-001,  -5.00000E-001,
 -1.00000E+000,   5.00000E-001,   5.00000E-001
};

const int pylith::meshio::MeshDataLagritTet::_cells[] = {
  9,      5,      8,     10,
  8,      9,      1,      5,
  1,      3,      2,      4,
  8,      5,      4,     10,
  8,      1,      4,      5,
  5,      3,      6,      7,
  4,      3,      2,      6,
  5,      3,      4,      6,
  1,      3,      4,      5,
  9,      5,     10,     11,
  8,      1,      0,      4,
  0,      1,      2,      4
};
const int pylith::meshio::MeshDataLagritTet::_materialIds[] = {
  2, 2, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1
};

const int pylith::meshio::MeshDataLagritTet::_numGroups = 7;

const int pylith::meshio::MeshDataLagritTet::_groupSizes[] = 
  { 4, 4, 4, 6, 6, 6, 6 };

const int pylith::meshio::MeshDataLagritTet::_groups[] = {
  0,  1,  4,  5,
  8,  9, 10, 11,
  2,  3,  6,  7,
  0,  1,  2,  3,  8,  9,
  4,  5,  6,  7, 10, 11,
  0,  2,  4,  6,  8, 10,
  1,  3,  5,  7,  9, 11
};

const char* pylith::meshio::MeshDataLagritTet::_groupNames[] = {
  "fault", "xm", "xp", "ym", "yp", "zm", "zp"
};

const char* pylith::meshio::MeshDataLagritTet::_groupTypes[] = {
  "vertex", "vertex", "vertex", "vertex", "vertex", "vertex", "vertex"
};

const bool pylith::meshio::MeshDataLagritTet::_useIndexZero = true;

pylith::meshio::MeshDataLagritTet::MeshDataLagritTet(void)
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

pylith::meshio::MeshDataLagritTet::~MeshDataLagritTet(void)
{}


// End of file
