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

#include "TestDataWriterHDF5ExtSubMeshCases.hh" // Implementation of class methods

#include "data/DataWriterHDF5DataSubMeshLine2.hh"
#include "data/DataWriterHDF5DataSubMeshTri3.hh"
#include "data/DataWriterHDF5DataSubMeshQuad4.hh"
#include "data/DataWriterHDF5DataSubMeshTet4.hh"
#include "data/DataWriterHDF5DataSubMeshHex8.hh"


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtSubMeshLine2 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtSubMeshTri3 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtSubMeshQuad4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtSubMeshTet4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtSubMeshHex8 );


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtSubMeshLine2::setUp(void)
{ // setUp
  TestDataWriterHDF5ExtSubMesh::setUp();
  _data = new DataWriterHDF5DataSubMeshLine2;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtSubMeshTri3::setUp(void)
{ // setUp
  TestDataWriterHDF5ExtSubMesh::setUp();
  _data = new DataWriterHDF5DataSubMeshTri3;
  _flipFault = true;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtSubMeshQuad4::setUp(void)
{ // setUp
  TestDataWriterHDF5ExtSubMesh::setUp();
  _data = new DataWriterHDF5DataSubMeshQuad4;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtSubMeshTet4::setUp(void)
{ // setUp
  TestDataWriterHDF5ExtSubMesh::setUp();
  _data = new DataWriterHDF5DataSubMeshTet4;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtSubMeshHex8::setUp(void)
{ // setUp
  TestDataWriterHDF5ExtSubMesh::setUp();
  _data = new DataWriterHDF5DataSubMeshHex8;
  _flipFault = true;
  _initialize();
} // setUp


// End of file 
