// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include "MeshData3D.hh"

const int pylith::meshio::MeshData3D::_numVertices = 14;

const int pylith::meshio::MeshData3D::_spaceDim = 3;

const int pylith::meshio::MeshData3D::_numCells = 2;

const int pylith::meshio::MeshData3D::_cellDim = 3;

const int pylith::meshio::MeshData3D::_numCorners = 8;

const PylithScalar pylith::meshio::MeshData3D::_vertices[] = {
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

const int pylith::meshio::MeshData3D::_numGroups = 2;

const int pylith::meshio::MeshData3D::_groupSizes[] = 
  { 5, 2, 4 };

const int pylith::meshio::MeshData3D::_groups[] = {
  0, 4, 6, 7, 10,
  0, 1,
  0, 4, 12, 13
};

const char* pylith::meshio::MeshData3D::_groupNames[] = {
  "group A", "group B", "group C"
};

const char* pylith::meshio::MeshData3D::_groupTypes[] = {
  "vertex", "cell", "vertex"
};

const bool pylith::meshio::MeshData3D::_useIndexZero = true;

pylith::meshio::MeshData3D::MeshData3D(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numCorners = _numCorners;
  vertices = const_cast<PylithScalar*>(_vertices);
  cells = const_cast<int*>(_cells);
  materialIds = const_cast<int*>(_materialIds);
  groups = const_cast<int*>(_groups);
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
  useIndexZero = _useIndexZero;
} // constructor

pylith::meshio::MeshData3D::~MeshData3D(void)
{}


// End of file
