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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterHDF5FaultMeshCases.hh" // Implementation of class methods

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END


// ----------------------------------------------------------------------
#include "data/DataWriterHDF5DataFaultMeshTri3.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5FaultMeshTri3 );

// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5FaultMeshTri3::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5FaultMesh::setUp();
  _data = new DataWriterHDF5DataFaultMeshTri3;

  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/DataWriterHDF5DataFaultMeshQuad4.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5FaultMeshQuad4 );

// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5FaultMeshQuad4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5FaultMesh::setUp();
  _data = new DataWriterHDF5DataFaultMeshQuad4;

  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/DataWriterHDF5DataFaultMeshTet4.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5FaultMeshTet4 );

// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5FaultMeshTet4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5FaultMesh::setUp();
  _data = new DataWriterHDF5DataFaultMeshTet4;

  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/DataWriterHDF5DataFaultMeshHex8.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5FaultMeshHex8 );

// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5FaultMeshHex8::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5FaultMesh::setUp();
  _data = new DataWriterHDF5DataFaultMeshHex8;

  _initialize();

  PYLITH_METHOD_END;
} // setUp


// End of file 
