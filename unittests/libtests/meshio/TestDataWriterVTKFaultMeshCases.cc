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

#include "TestDataWriterVTKFaultMeshCases.hh" // Implementation of class methods

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
#include "data/DataWriterVTKDataFaultMeshTri3.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKFaultMeshTri3 );

// Setup testing data.
void
pylith::meshio::TestDataWriterVTKFaultMeshTri3::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterVTKFaultMesh::setUp();
  _data = new DataWriterVTKDataFaultMeshTri3;

  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/DataWriterVTKDataFaultMeshQuad4.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKFaultMeshQuad4 );

// Setup testing data.
void
pylith::meshio::TestDataWriterVTKFaultMeshQuad4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterVTKFaultMesh::setUp();
  _data = new DataWriterVTKDataFaultMeshQuad4;

  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/DataWriterVTKDataFaultMeshTet4.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKFaultMeshTet4 );

// Setup testing data.
void
pylith::meshio::TestDataWriterVTKFaultMeshTet4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterVTKFaultMesh::setUp();
  _data = new DataWriterVTKDataFaultMeshTet4;

  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/DataWriterVTKDataFaultMeshHex8.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKFaultMeshHex8 );

// Setup testing data.
void
pylith::meshio::TestDataWriterVTKFaultMeshHex8::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterVTKFaultMesh::setUp();
  _data = new DataWriterVTKDataFaultMeshHex8;

  _initialize();

  PYLITH_METHOD_END;
} // setUp


// End of file 
