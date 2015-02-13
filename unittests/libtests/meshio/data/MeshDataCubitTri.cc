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

#include "MeshDataCubitTri.hh"

const int pylith::meshio::MeshDataCubitTri::_numVertices = 4;

const int pylith::meshio::MeshDataCubitTri::_spaceDim = 2;

const int pylith::meshio::MeshDataCubitTri::_numCells = 2;

const int pylith::meshio::MeshDataCubitTri::_cellDim = 2;

const int pylith::meshio::MeshDataCubitTri::_numCorners = 3;

const PylithScalar pylith::meshio::MeshDataCubitTri::_vertices[] = {
  -1.0,  0.0,
   0.0, -1.0,
   0.0,  1.0,
   1.0,  0.0
};

const int pylith::meshio::MeshDataCubitTri::_cells[] = {
  0,  1,  2,
  2,  1,  3
};
const int pylith::meshio::MeshDataCubitTri::_materialIds[] = {
  2, 3
};

const int pylith::meshio::MeshDataCubitTri::_numGroups = 2;

const int pylith::meshio::MeshDataCubitTri::_groupSizes[] = 
  { 1, 2 };

const int pylith::meshio::MeshDataCubitTri::_groups[] = {
  0,
  2, 3
};

const char* pylith::meshio::MeshDataCubitTri::_groupNames[] = {
  "left_vertex", "right_vertex"
};

const char* pylith::meshio::MeshDataCubitTri::_groupTypes[] = {
  "vertex", "vertex"
};

const bool pylith::meshio::MeshDataCubitTri::_useIndexZero = true;

pylith::meshio::MeshDataCubitTri::MeshDataCubitTri(void)
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

pylith::meshio::MeshDataCubitTri::~MeshDataCubitTri(void)
{}


// End of file
