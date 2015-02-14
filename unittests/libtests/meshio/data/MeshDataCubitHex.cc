// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "MeshDataCubitHex.hh"

const int pylith::meshio::MeshDataCubitHex::_numVertices = 12;

const int pylith::meshio::MeshDataCubitHex::_spaceDim = 3;

const int pylith::meshio::MeshDataCubitHex::_numCells = 2;

const int pylith::meshio::MeshDataCubitHex::_cellDim = 3;

const int pylith::meshio::MeshDataCubitHex::_numCorners = 8;

const PylithScalar pylith::meshio::MeshDataCubitHex::_vertices[] = {
  -2.0, -1.0,  1.0,
  -2.0, -1.0, -1.0,
  -2.0,  1.0, -1.0,
  -2.0,  1.0,  1.0,
   0.0, -1.0,  1.0,
   0.0, -1.0, -1.0,
   0.0,  1.0, -1.0,
   0.0,  1.0,  1.0,
   2.0, -1.0,  1.0,
   2.0, -1.0, -1.0,
   2.0,  1.0, -1.0,
   2.0,  1.0,  1.0
};

// PyLith order, not Cubit order
const int pylith::meshio::MeshDataCubitHex::_cells[] = {
  0,  1,  2,  3,  4,  5,  6,  7,
  4,  5,  6,  7,  8,  9, 10, 11
};
const int pylith::meshio::MeshDataCubitHex::_materialIds[] = {
  7, 8
};

const int pylith::meshio::MeshDataCubitHex::_numGroups = 2;

const int pylith::meshio::MeshDataCubitHex::_groupSizes[] = 
  { 4, 6 };

const int pylith::meshio::MeshDataCubitHex::_groups[] = {
  8,  9, 10, 11,
  0,  3,  4,  7,  8, 11
};

const char* pylith::meshio::MeshDataCubitHex::_groupNames[] = {
  "right_face", "top_face"
};

const char* pylith::meshio::MeshDataCubitHex::_groupTypes[] = {
  "vertex", "vertex"
};

const bool pylith::meshio::MeshDataCubitHex::_useIndexZero = true;

pylith::meshio::MeshDataCubitHex::MeshDataCubitHex(void)
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

pylith::meshio::MeshDataCubitHex::~MeshDataCubitHex(void)
{}


// End of file
