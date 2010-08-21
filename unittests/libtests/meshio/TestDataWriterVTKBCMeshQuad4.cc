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

#include "TestDataWriterVTKBCMeshQuad4.hh" // Implementation of class methods

#include "data/DataWriterVTKDataBCMeshQuad4.hh" // USES DataWriterVTKDataMeshQuad4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKBCMeshQuad4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKBCMeshQuad4::setUp(void)
{ // setUp
  TestDataWriterVTKBCMesh::setUp();
  _data = new DataWriterVTKDataBCMeshQuad4;
  _initialize();
} // setUp


// End of file 
