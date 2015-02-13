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

#include "BoundaryMeshDataHex8.hh"

const char* pylith::bc::BoundaryMeshDataHex8::_filename = "data/hex8.mesh";

const char* pylith::bc::BoundaryMeshDataHex8::_bcLabel = "bc2";

const char* pylith::bc::BoundaryMeshDataHex8::_faultLabel = "fault";
const int pylith::bc::BoundaryMeshDataHex8::_faultId = 100;

const int pylith::bc::BoundaryMeshDataHex8::_numCorners = 4;
const int pylith::bc::BoundaryMeshDataHex8::_numCells = 2;
const bool pylith::bc::BoundaryMeshDataHex8::_isSimplexMesh = false;

const int pylith::bc::BoundaryMeshDataHex8::_numVerticesNoFault = 6;
const int pylith::bc::BoundaryMeshDataHex8::_numVerticesFault = 8;

pylith::bc::BoundaryMeshDataHex8::BoundaryMeshDataHex8(void)
{ // constructor
  filename = const_cast<char*>(_filename);

  bcLabel = const_cast<char*>(_bcLabel);

  faultLabel = const_cast<char*>(_faultLabel);
  faultId = _faultId;

  numCorners = _numCorners;
  numCells = _numCells;
  isSimplexMesh = _isSimplexMesh;

  numVerticesNoFault = _numVerticesNoFault;

  numVerticesFault = _numVerticesFault;
} // constructor

pylith::bc::BoundaryMeshDataHex8::~BoundaryMeshDataHex8(void)
{}


// End of file
