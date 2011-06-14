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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterHDF5ExtFaultMeshCases.hh" // Implementation of class methods

#include "data/DataWriterHDF5DataFaultMeshTri3.hh"
#include "data/DataWriterHDF5DataFaultMeshQuad4.hh"
#include "data/DataWriterHDF5DataFaultMeshTet4.hh"
#include "data/DataWriterHDF5DataFaultMeshHex8.hh"


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtFaultMeshTri3 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtFaultMeshQuad4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtFaultMeshTet4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtFaultMeshHex8 );


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtFaultMeshTri3::setUp(void)
{ // setUp
  TestDataWriterHDF5ExtFaultMesh::setUp();
  _data = new DataWriterHDF5DataFaultMeshTri3;
  _flipFault = true;

  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtFaultMeshQuad4::setUp(void)
{ // setUp
  TestDataWriterHDF5ExtFaultMesh::setUp();
  _data = new DataWriterHDF5DataFaultMeshQuad4;
  _flipFault = true;

  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtFaultMeshTet4::setUp(void)
{ // setUp
  TestDataWriterHDF5ExtFaultMesh::setUp();
  _data = new DataWriterHDF5DataFaultMeshTet4;
  _flipFault = false;

  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtFaultMeshHex8::setUp(void)
{ // setUp
  TestDataWriterHDF5ExtFaultMesh::setUp();
  _data = new DataWriterHDF5DataFaultMeshHex8;
  _flipFault = true;

  _initialize();
} // setUp


// End of file 
