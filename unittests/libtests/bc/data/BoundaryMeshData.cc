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

#include "BoundaryMeshData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::BoundaryMeshData::BoundaryMeshData(void) :
  filename(0),
  bcLabel(0),
  faultLabel(0),
  faultId(0),
  isSimplexMesh(true)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::bc::BoundaryMeshData::~BoundaryMeshData(void)
{ // destructor
} // destructor


// End of file
