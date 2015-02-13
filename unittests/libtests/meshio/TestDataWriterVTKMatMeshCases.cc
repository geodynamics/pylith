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

#include "TestDataWriterVTKMatMeshCases.hh" // Implementation of class methods

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include "data/DataWriterVTKDataMatMeshTri3.hh"
#include "data/DataWriterVTKDataMatMeshQuad4.hh"
#include "data/DataWriterVTKDataMatMeshTet4.hh"
#include "data/DataWriterVTKDataMatMeshHex8.hh"


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKMatMeshTri3 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKMatMeshQuad4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKMatMeshTet4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKMatMeshHex8 );


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKMatMeshTri3::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterVTKMesh::setUp();
  _data = new DataWriterVTKDataMatMeshTri3;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKMatMeshQuad4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterVTKMesh::setUp();
  _data = new DataWriterVTKDataMatMeshQuad4;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKMatMeshTet4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterVTKMesh::setUp();
  _data = new DataWriterVTKDataMatMeshTet4;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKMatMeshHex8::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterVTKMesh::setUp();
  _data = new DataWriterVTKDataMatMeshHex8;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// End of file 
