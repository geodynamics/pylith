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

#include "TestDataWriterHDF5ExtMeshHex8.hh" // Implementation of class methods

#include "data/DataWriterHDF5DataMeshHex8.hh" // USES DataWriterHDF5DataMeshHex8

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtMeshHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMeshHex8::setUp(void)
{ // setUp
  TestDataWriterHDF5ExtMesh::setUp();
  _data = new DataWriterHDF5DataMeshHex8;
  _flipFault = true;

  _initialize();
} // setUp


// End of file 
