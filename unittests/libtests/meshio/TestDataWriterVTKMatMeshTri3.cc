// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterVTKMatMeshTri3.hh" // Implementation of class methods

#include "data/DataWriterVTKDataMatMeshTri3.hh" // USES DataWriterVTKDataMeshTri3

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKMatMeshTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKMatMeshTri3::setUp(void)
{ // setUp
  TestDataWriterVTKMesh::setUp();
  _data = new DataWriterVTKDataMatMeshTri3;
  _flipFault = true;

  _initialize();
} // setUp


// End of file 
