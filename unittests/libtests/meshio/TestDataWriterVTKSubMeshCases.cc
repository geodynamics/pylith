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

#include "TestDataWriterVTKSubMeshCases.hh" // Implementation of class methods

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END


// ----------------------------------------------------------------------
#include "data/DataWriterVTKDataSubMeshTri3.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKSubMeshTri3 );

// Setup testing data.
void
pylith::meshio::TestDataWriterVTKSubMeshTri3::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterVTKSubMesh::setUp();
  _data = new DataWriterVTKDataSubMeshTri3;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/DataWriterVTKDataSubMeshQuad4.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKSubMeshQuad4 );

// Setup testing data.
void
pylith::meshio::TestDataWriterVTKSubMeshQuad4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterVTKSubMesh::setUp();
  _data = new DataWriterVTKDataSubMeshQuad4;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/DataWriterVTKDataSubMeshTet4.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKSubMeshTet4 );

// Setup testing data.
void
pylith::meshio::TestDataWriterVTKSubMeshTet4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterVTKSubMesh::setUp();
  _data = new DataWriterVTKDataSubMeshTet4;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/DataWriterVTKDataSubMeshHex8.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKSubMeshHex8 );

// Setup testing data.
void
pylith::meshio::TestDataWriterVTKSubMeshHex8::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterVTKSubMesh::setUp();
  _data = new DataWriterVTKDataSubMeshHex8;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// End of file 
