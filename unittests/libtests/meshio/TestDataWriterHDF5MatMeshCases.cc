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

#include "TestDataWriterHDF5MatMeshCases.hh" // Implementation of class methods

#include "data/DataWriterHDF5DataMatMeshTri3.hh"
#include "data/DataWriterHDF5DataMatMeshQuad4.hh"
#include "data/DataWriterHDF5DataMatMeshTet4.hh"
#include "data/DataWriterHDF5DataMatMeshHex8.hh"

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5MatMeshTri3 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5MatMeshQuad4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5MatMeshTet4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5MatMeshHex8 );


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5MatMeshTri3::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5Mesh::setUp();
  _data = new DataWriterHDF5DataMatMeshTri3;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5MatMeshQuad4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5Mesh::setUp();
  _data = new DataWriterHDF5DataMatMeshQuad4;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5MatMeshTet4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5Mesh::setUp();
  _data = new DataWriterHDF5DataMatMeshTet4;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5MatMeshHex8::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5Mesh::setUp();
  _data = new DataWriterHDF5DataMatMeshHex8;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// End of file 
