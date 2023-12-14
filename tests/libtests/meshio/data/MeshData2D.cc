// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
nst int pylith::meshio::MeshData2D::_numCells = 3;

const int pylith::meshio::MeshData2D::_cellDim = 2;

const int pylith::meshio::MeshData2D::_numCorners = 4;

const PylithScalar pylith::meshio::MeshData2D::_vertices[] = {
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
const int pylith::meshio::MeshData2D::_materialIds[] = {
    1, 0, 1
};

const int pylith::meshio::MeshData2D::_numGroups = 3;

const int pylith::meshio::MeshData2D::_groupSizes[] =
{ 5, 3, 2 };

const int pylith::meshio::MeshData2D::_groups[] = {
    0, 2, 4, 6, 8,
    1, 4, 7,
    0, 2
};

const char* pylith::meshio::MeshData2D::_groupNames[] = {
    "group A", "group B", "group C"
};

const char* pylith::meshio::MeshData2D::_groupTypes[] = {
    "vertex", "vertex", "cell"
};

const bool pylith::meshio::MeshData2D::_useIndexZero = true;

pylith::meshio::MeshData2D::MeshData2D(void) { // constructor
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


pylith::meshio::MeshData2D::~MeshData2D(void) {}


// End of file
