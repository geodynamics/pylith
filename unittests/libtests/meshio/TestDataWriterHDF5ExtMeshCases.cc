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

#include "TestDataWriterHDF5ExtMeshCases.hh" // Implementation of class methods

#include "data/DataWriterHDF5DataMeshTri3.hh"
#include "data/DataWriterHDF5DataMeshQuad4.hh"
#include "data/DataWriterHDF5DataMeshTet4.hh"
#include "data/DataWriterHDF5DataMeshHex8.hh"

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtMeshTri3 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtMeshQuad4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtMeshTet4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtMeshHex8 );


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMeshTri3::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5ExtMesh::setUp();
  _data = new DataWriterHDF5DataMeshTri3;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMeshQuad4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5ExtMesh::setUp();
  _data = new DataWriterHDF5DataMeshQuad4;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMeshTet4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5ExtMesh::setUp();
  _data = new DataWriterHDF5DataMeshTet4;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMeshHex8::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5ExtMesh::setUp();
  _data = new DataWriterHDF5DataMeshHex8;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// End of file 
