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

#include "MeshData1Din2D.hh"

const int pylith::meshio::MeshData1Din2D::_numVertices = 4;

const int pylith::meshio::MeshData1Din2D::_spaceDim = 2;

const int pylith::meshio::MeshData1Din2D::_numCells = 3;

const int pylith::meshio::MeshData1Din2D::_cellDim = 1;

const int pylith::meshio::MeshData1Din2D::_numCorners = 2;

const PylithScalar pylith::meshio::MeshData1Din2D::_vertices[] = {
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

const int pylith::meshio::MeshData1Din2D::_numGroups = 2;

const int pylith::meshio::MeshData1Din2D::_groupSizes[] = 
  { 2, 3 };

const int pylith::meshio::MeshData1Din2D::_groups[] = {
  0, 2,
  0, 1, 3
};

const char* pylith::meshio::MeshData1Din2D::_groupNames[] = {
  "group A", "group B"
};

const char* pylith::meshio::MeshData1Din2D::_groupTypes[] = {
  "cell", "vertex"
};

const bool pylith::meshio::MeshData1Din2D::_useIndexZero = true;

pylith::meshio::MeshData1Din2D::MeshData1Din2D(void)
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

pylith::meshio::MeshData1Din2D::~MeshData1Din2D(void)
{}


// End of file
