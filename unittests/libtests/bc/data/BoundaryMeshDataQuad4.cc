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
// Copyright (c) 2010-2013 University of California, Davis
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

const int pylith::bc::BoundaryMeshDataQuad4::_numVerticesNoFault = 3;
const int pylith::bc::BoundaryMeshDataQuad4::_verticesNoFault[] = {
  2, 4, 6
};
const int pylith::bc::BoundaryMeshDataQuad4::_cellsNoFault[] = {
  2, 4,
  4, 6,
};

const int pylith::bc::BoundaryMeshDataQuad4::_numVerticesFault = 4;
const int pylith::bc::BoundaryMeshDataQuad4::_verticesFault[] = {
  3, 5, 7, 9,
};
const int pylith::bc::BoundaryMeshDataQuad4::_cellsFault[] = {
  3, 9,
  5, 7,
};

pylith::bc::BoundaryMeshDataQuad4::BoundaryMeshDataQuad4(void)
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

pylith::bc::BoundaryMeshDataQuad4::~BoundaryMeshDataQuad4(void)
{}


// End of file
