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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterHDF5SubMeshCases.hh" // Implementation of class methods

#include "data/DataWriterHDF5DataSubMeshLine2.hh"
#include "data/DataWriterHDF5DataSubMeshTri3.hh"
#include "data/DataWriterHDF5DataSubMeshQuad4.hh"
#include "data/DataWriterHDF5DataSubMeshTet4.hh"
#include "data/DataWriterHDF5DataSubMeshHex8.hh"


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5SubMeshLine2 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5SubMeshTri3 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5SubMeshQuad4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5SubMeshTet4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5SubMeshHex8 );


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5SubMeshLine2::setUp(void)
{ // setUp
  TestDataWriterHDF5SubMesh::setUp();
  _data = new DataWriterHDF5DataSubMeshLine2;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5SubMeshTri3::setUp(void)
{ // setUp
  TestDataWriterHDF5SubMesh::setUp();
  _data = new DataWriterHDF5DataSubMeshTri3;
  _flipFault = true;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5SubMeshQuad4::setUp(void)
{ // setUp
  TestDataWriterHDF5SubMesh::setUp();
  _data = new DataWriterHDF5DataSubMeshQuad4;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5SubMeshTet4::setUp(void)
{ // setUp
  TestDataWriterHDF5SubMesh::setUp();
  _data = new DataWriterHDF5DataSubMeshTet4;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5SubMeshHex8::setUp(void)
{ // setUp
  TestDataWriterHDF5SubMesh::setUp();
  _data = new DataWriterHDF5DataSubMeshHex8;
  _flipFault = true;
  _initialize();
} // setUp


// End of file 
