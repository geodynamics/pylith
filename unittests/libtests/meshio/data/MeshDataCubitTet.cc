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

#include "MeshDataCubitTet.hh"

const int pylith::meshio::MeshDataCubitTet::_numVertices = 5;

const int pylith::meshio::MeshDataCubitTet::_spaceDim = 3;

const int pylith::meshio::MeshDataCubitTet::_numCells = 2;

const int pylith::meshio::MeshDataCubitTet::_cellDim = 3;

const int pylith::meshio::MeshDataCubitTet::_numCorners = 4;

const PylithScalar pylith::meshio::MeshDataCubitTet::_vertices[] = {
  -2.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  1.0,  0.0,
   0.0,  0.0,  2.0,
   2.0,  0.0,  0.0
};

const int pylith::meshio::MeshDataCubitTet::_cells[] = {
  0,  1,  2,  3,
  1,  4,  2,  3
};
const int pylith::meshio::MeshDataCubitTet::_materialIds[] = {
  7, 8
};

const int pylith::meshio::MeshDataCubitTet::_numGroups = 2;

const int pylith::meshio::MeshDataCubitTet::_groupSizes[] = 
  { 3, 4 };

const int pylith::meshio::MeshDataCubitTet::_groups[] = {
  1, 2, 3,
  0, 1, 2, 3,
};

const char* pylith::meshio::MeshDataCubitTet::_groupNames[] = {
  "mid_face",
  "bottom_face",
};

const char* pylith::meshio::MeshDataCubitTet::_groupTypes[] = {
  "vertex", "vertex"
};

const bool pylith::meshio::MeshDataCubitTet::_useIndexZero = true;

pylith::meshio::MeshDataCubitTet::MeshDataCubitTet(void)
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

pylith::meshio::MeshDataCubitTet::~MeshDataCubitTet(void)
{}


// End of file
