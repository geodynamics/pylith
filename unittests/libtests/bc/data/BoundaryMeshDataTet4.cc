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

#include "BoundaryMeshDataTet4.hh"

const char* pylith::bc::BoundaryMeshDataTet4::_filename = "data/tet4.mesh";

const char* pylith::bc::BoundaryMeshDataTet4::_bcLabel = "bc4";

const char* pylith::bc::BoundaryMeshDataTet4::_faultLabel = "fault";
const int pylith::bc::BoundaryMeshDataTet4::_faultId = 100;

const int pylith::bc::BoundaryMeshDataTet4::_numCorners = 3;
const int pylith::bc::BoundaryMeshDataTet4::_numCells = 2;

const int pylith::bc::BoundaryMeshDataTet4::_numVerticesNoFault = 4;
const int pylith::bc::BoundaryMeshDataTet4::_verticesNoFault[] = {
  2, 4, 5, 6
};
const int pylith::bc::BoundaryMeshDataTet4::_cellsNoFault[] = {
  5, 4, 2,
  5, 2, 6,
};

const int pylith::bc::BoundaryMeshDataTet4::_numVerticesFault = 6;
const int pylith::bc::BoundaryMeshDataTet4::_verticesFault[] = {
  3, 5, 6, 7, 8, 10
};
const int pylith::bc::BoundaryMeshDataTet4::_cellsFault[] = {
  6, 5, 3,
  10, 8, 7,
};

pylith::bc::BoundaryMeshDataTet4::BoundaryMeshDataTet4(void)
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

pylith::bc::BoundaryMeshDataTet4::~BoundaryMeshDataTet4(void)
{}


// End of file
