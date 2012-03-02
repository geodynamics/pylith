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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterVTKPointsCases.hh" // Implementation of class methods

#include "data/DataWriterVTKDataPointsTri3.hh" // USES DataWriterVTKDataPointsTri3
#include "data/DataWriterVTKDataPointsQuad4.hh" // USES DataWriterVTKDataPointsQuad4
#include "data/DataWriterVTKDataPointsTet4.hh" // USES DataWriterVTKDataPointsTet4
#include "data/DataWriterVTKDataPointsHex8.hh" // USES DataWriterVTKDataPointsHex8


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKPointsTri3 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKPointsQuad4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKPointsTet4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKPointsHex8 );


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKPointsTri3::setUp(void)
{ // setUp
  TestDataWriterVTKPoints::setUp();
  _data = new DataWriterVTKDataPointsTri3;
  _flipFault = true;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKPointsQuad4::setUp(void)
{ // setUp
  TestDataWriterVTKPoints::setUp();
  _data = new DataWriterVTKDataPointsQuad4;
  _flipFault = false;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKPointsTet4::setUp(void)
{ // setUp
  TestDataWriterVTKPoints::setUp();
  _data = new DataWriterVTKDataPointsTet4;
  _flipFault = true;
  _initialize();
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKPointsHex8::setUp(void)
{ // setUp
  TestDataWriterVTKPoints::setUp();
  _data = new DataWriterVTKDataPointsHex8;
  _flipFault = true;
  _initialize();
} // setUp


// End of file 
