// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include "MeshData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshData::MeshData(void) :
  numVertices(0),
  spaceDim(0),
  numCells(0),
  cellDim(0),
  numCorners(0),
  vertices(0),
  cells(0),
  materialIds(0),
  groups(0),
  groupSizes(0),
  groupNames(0),
  groupTypes(0),
  numGroups(0),
  useIndexZero(true)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshData::~MeshData(void)
{ // destructor
} // destructor


// End of file
