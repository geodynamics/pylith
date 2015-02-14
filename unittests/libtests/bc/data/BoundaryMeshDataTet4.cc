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

#include "BoundaryMeshDataTet4.hh"

const char* pylith::bc::BoundaryMeshDataTet4::_filename = "data/tet4.mesh";

const char* pylith::bc::BoundaryMeshDataTet4::_bcLabel = "bc4";

const char* pylith::bc::BoundaryMeshDataTet4::_faultLabel = "fault";
const int pylith::bc::BoundaryMeshDataTet4::_faultId = 100;

const int pylith::bc::BoundaryMeshDataTet4::_numCorners = 3;
const int pylith::bc::BoundaryMeshDataTet4::_numCells = 2;
const bool pylith::bc::BoundaryMeshDataTet4::_isSimplexMesh = false;

const int pylith::bc::BoundaryMeshDataTet4::_numVerticesNoFault = 4;

const int pylith::bc::BoundaryMeshDataTet4::_numVerticesFault = 6;

pylith::bc::BoundaryMeshDataTet4::BoundaryMeshDataTet4(void)
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

pylith::bc::BoundaryMeshDataTet4::~BoundaryMeshDataTet4(void)
{}


// End of file
