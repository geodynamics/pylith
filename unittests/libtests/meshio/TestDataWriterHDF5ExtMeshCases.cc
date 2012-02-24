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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterHDF5ExtMeshCases.hh" // Implementation of class methods

#include "data/DataWriterHDF5DataMeshLine2.hh"
#include "data/DataWriterHDF5DataMeshTri3.hh"
#include "data/DataWriterHDF5DataMeshQuad4.hh"
#include "data/DataWriterHDF5DataMeshTet4.hh"
#include "data/DataWriterHDF5DataMeshHex8.hh"


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtMeshLine2 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtMeshTri3 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtMeshQuad4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtMeshTet4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtMeshHex8 );


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMeshLine2::setUp(void)
{ // setUp
  TestDataWriterHDF5ExtMesh::setUp();
  _data = new DataWriterHDF5DataMeshLine2;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMeshTri3::setUp(void)
{ // setUp
  TestDataWriterHDF5ExtMesh::setUp();
  _data = new DataWriterHDF5DataMeshTri3;
  _flipFault = true;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMeshQuad4::setUp(void)
{ // setUp
  TestDataWriterHDF5ExtMesh::setUp();
  _data = new DataWriterHDF5DataMeshQuad4;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMeshTet4::setUp(void)
{ // setUp
  TestDataWriterHDF5ExtMesh::setUp();
  _data = new DataWriterHDF5DataMeshTet4;
  _flipFault = false;
  _initialize();
} // setUp


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
