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

#include "BoundaryMeshDataTri3.hh"

const char* pylith::bc::BoundaryMeshDataTri3::_filename = "data/tri3.mesh";

const char* pylith::bc::BoundaryMeshDataTri3::_bcLabel = "bc";

const char* pylith::bc::BoundaryMeshDataTri3::_faultLabel = "fault";
const int pylith::bc::BoundaryMeshDataTri3::_faultId = 100;

const int pylith::bc::BoundaryMeshDataTri3::_numCorners = 2;
const int pylith::bc::BoundaryMeshDataTri3::_numCells = 1;
const bool pylith::bc::BoundaryMeshDataTri3::_isSimplexMesh = false;

const int pylith::bc::BoundaryMeshDataTri3::_numVerticesNoFault = 2;

const int pylith::bc::BoundaryMeshDataTri3::_numVerticesFault = 2;

pylith::bc::BoundaryMeshDataTri3::BoundaryMeshDataTri3(void)
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

pylith::bc::BoundaryMeshDataTri3::~BoundaryMeshDataTri3(void)
{}


// End of file
