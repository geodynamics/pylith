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

#include "TestDataWriterVTKFaultMeshHex8.hh" // Implementation of class methods

#include "data/DataWriterVTKDataFaultMeshHex8.hh" // USES DataWriterVTKDataFaultMeshHex8

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKFaultMeshHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKFaultMeshHex8::setUp(void)
{ // setUp
  TestDataWriterVTKFaultMesh::setUp();
  _data = new DataWriterVTKDataFaultMeshHex8;
  _flipFault = true;

  _initialize();
} // setUp


// End of file 
