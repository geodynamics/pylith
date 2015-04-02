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

#include "TestDataWriterHDF5SubMeshCases.hh" // Implementation of class methods

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END


// ----------------------------------------------------------------------
#include "data/DataWriterHDF5DataSubMeshTri3.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5SubMeshTri3 );

// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5SubMeshTri3::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5SubMesh::setUp();
  _data = new DataWriterHDF5DataSubMeshTri3;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/DataWriterHDF5DataSubMeshQuad4.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5SubMeshQuad4 );

// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5SubMeshQuad4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5SubMesh::setUp();
  _data = new DataWriterHDF5DataSubMeshQuad4;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/DataWriterHDF5DataSubMeshTet4.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5SubMeshTet4 );

// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5SubMeshTet4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5SubMesh::setUp();
  _data = new DataWriterHDF5DataSubMeshTet4;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/DataWriterHDF5DataSubMeshHex8.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5SubMeshHex8 );

// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5SubMeshHex8::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5SubMesh::setUp();
  _data = new DataWriterHDF5DataSubMeshHex8;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// End of file 
