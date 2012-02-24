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

#include "TestDataWriterHDF5MatMeshCases.hh" // Implementation of class methods

#include "data/DataWriterHDF5DataMatMeshLine2.hh"
#include "data/DataWriterHDF5DataMatMeshTri3.hh"
#include "data/DataWriterHDF5DataMatMeshQuad4.hh"
#include "data/DataWriterHDF5DataMatMeshTet4.hh"
#include "data/DataWriterHDF5DataMatMeshHex8.hh"


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5MatMeshLine2 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5MatMeshTri3 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5MatMeshQuad4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5MatMeshTet4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5MatMeshHex8 );


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5MatMeshLine2::setUp(void)
{ // setUp
  TestDataWriterHDF5Mesh::setUp();
  _data = new DataWriterHDF5DataMatMeshLine2;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5MatMeshTri3::setUp(void)
{ // setUp
  TestDataWriterHDF5Mesh::setUp();
  _data = new DataWriterHDF5DataMatMeshTri3;
  _flipFault = true;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5MatMeshQuad4::setUp(void)
{ // setUp
  TestDataWriterHDF5Mesh::setUp();
  _data = new DataWriterHDF5DataMatMeshQuad4;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5MatMeshTet4::setUp(void)
{ // setUp
  TestDataWriterHDF5Mesh::setUp();
  _data = new DataWriterHDF5DataMatMeshTet4;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5MatMeshHex8::setUp(void)
{ // setUp
  TestDataWriterHDF5Mesh::setUp();
  _data = new DataWriterHDF5DataMatMeshHex8;
  _flipFault = true;
  _initialize();
} // setUp


// End of file 
