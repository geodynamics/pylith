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

#include "MeshDataCohesive.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::topology::MeshDataCohesive::MeshDataCohesive(void) :
  filename(0),
  refineLevel(0),
  faultA(0),
  faultB(0),
  numVertices(0),
  spaceDim(0),
  cellDim(0),
  numCells(0),
  numCorners(0),
  numCellsCohesive(0),
  numCornersCohesive(0),
  matIdSum(0),
  groupSizes(0),
  groupNames(0),
  groupTypes(0),
  numGroups(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::topology::MeshDataCohesive::~MeshDataCohesive(void)
{ // destructor
} // destructor


// End of file
