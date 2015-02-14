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

#include "BoundaryMeshDataQuad4.hh"

const char* pylith::bc::BoundaryMeshDataQuad4::_filename = "data/quad4.mesh";

const char* pylith::bc::BoundaryMeshDataQuad4::_bcLabel = "bc3";

const char* pylith::bc::BoundaryMeshDataQuad4::_faultLabel = "fault";
const int pylith::bc::BoundaryMeshDataQuad4::_faultId = 100;

const int pylith::bc::BoundaryMeshDataQuad4::_numCorners = 2;
const int pylith::bc::BoundaryMeshDataQuad4::_numCells = 2;
const bool pylith::bc::BoundaryMeshDataQuad4::_isSimplexMesh = false;

const int pylith::bc::BoundaryMeshDataQuad4::_numVerticesNoFault = 3;

const int pylith::bc::BoundaryMeshDataQuad4::_numVerticesFault = 4;

pylith::bc::BoundaryMeshDataQuad4::BoundaryMeshDataQuad4(void)
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

pylith::bc::BoundaryMeshDataQuad4::~BoundaryMeshDataQuad4(void)
{}


// End of file
