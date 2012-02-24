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

#include "TestDataWriterHDF5MeshCases.hh" // Implementation of class methods

#include "data/DataWriterHDF5DataMeshLine2.hh"
#include "data/DataWriterHDF5DataMeshTri3.hh"
#include "data/DataWriterHDF5DataMeshQuad4.hh"
#include "data/DataWriterHDF5DataMeshTet4.hh"
#include "data/DataWriterHDF5DataMeshHex8.hh"


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5MeshLine2 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5MeshTri3 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5MeshQuad4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5MeshTet4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5MeshHex8 );


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5MeshLine2::setUp(void)
{ // setUp
  TestDataWriterHDF5Mesh::setUp();
  _data = new DataWriterHDF5DataMeshLine2;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5MeshTri3::setUp(void)
{ // setUp
  TestDataWriterHDF5Mesh::setUp();
  _data = new DataWriterHDF5DataMeshTri3;
  _flipFault = true;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5MeshQuad4::setUp(void)
{ // setUp
  TestDataWriterHDF5Mesh::setUp();
  _data = new DataWriterHDF5DataMeshQuad4;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5MeshTet4::setUp(void)
{ // setUp
  TestDataWriterHDF5Mesh::setUp();
  _data = new DataWriterHDF5DataMeshTet4;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5MeshHex8::setUp(void)
{ // setUp
  TestDataWriterHDF5Mesh::setUp();
  _data = new DataWriterHDF5DataMeshHex8;
  _flipFault = true;
  _initialize();
} // setUp


// End of file 
