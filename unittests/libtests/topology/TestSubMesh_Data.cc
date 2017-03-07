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

#include "TestSubMesh_Data.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::topology::TestSubMesh_Data::TestSubMesh_Data(void) :
    cellDim(0),
    numVertices(0),
    numCells(0),
    numCorners(0),
    cells(NULL),
    coordinates(NULL),
    label(NULL),
    groupSize(0),
    groupVertices(NULL),
    submeshNumCorners(0),
    submeshNumVertices(0),
    submeshVertices(NULL),
    submeshNumCells(0),
    submeshCells(NULL)
{ // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::topology::TestSubMesh_Data::~TestSubMesh_Data(void)
{ // destructor
} // destructor

// End of file
