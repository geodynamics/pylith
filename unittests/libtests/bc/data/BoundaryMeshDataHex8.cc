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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "BoundaryMeshDataHex8.hh"

const char* pylith::bc::BoundaryMeshDataHex8::_filename = "data/hex8.mesh";

const char* pylith::bc::BoundaryMeshDataHex8::_bcLabel = "bc2";

const char* pylith::bc::BoundaryMeshDataHex8::_faultLabel = "fault";
const int pylith::bc::BoundaryMeshDataHex8::_faultId = 100;

const int pylith::bc::BoundaryMeshDataHex8::_numCorners = 4;
const int pylith::bc::BoundaryMeshDataHex8::_numCells = 2;

const int pylith::bc::BoundaryMeshDataHex8::_numVerticesNoFault = 6;
const int pylith::bc::BoundaryMeshDataHex8::_verticesNoFault[] = {
  2, 4, 6, 8, 10, 12
};
const int pylith::bc::BoundaryMeshDataHex8::_cellsNoFault[] = {
  2, 4, 10, 8,
  4, 6, 12, 10,
};

const int pylith::bc::BoundaryMeshDataHex8::_numVerticesFault = 8;
const int pylith::bc::BoundaryMeshDataHex8::_verticesFault[] = {
  3, 5, 7, 9, 11, 13, 15, 17
};
const int pylith::bc::BoundaryMeshDataHex8::_cellsFault[] = {
  3, 15, 17, 9,
  5, 7, 13, 11,
};

pylith::bc::BoundaryMeshDataHex8::BoundaryMeshDataHex8(void)
{ // constructor
  filename = const_cast<char*>(_filename);

  bcLabel = const_cast<char*>(_bcLabel);

  faultLabel = const_cast<char*>(_faultLabel);
  faultId = _faultId;

  numCorners = _numCorners;
  numCells = _numCells;

  numVerticesNoFault = _numVerticesNoFault;
  verticesNoFault = const_cast<int*>(_verticesNoFault);
  cellsNoFault = const_cast<int*>(_cellsNoFault);

  numVerticesFault = _numVerticesFault;
  verticesFault = const_cast<int*>(_verticesFault);
  cellsFault = const_cast<int*>(_cellsFault);
} // constructor

pylith::bc::BoundaryMeshDataHex8::~BoundaryMeshDataHex8(void)
{}


// End of file
