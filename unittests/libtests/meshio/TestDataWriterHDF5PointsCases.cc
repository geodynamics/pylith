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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterHDF5PointsCases.hh" // Implementation of class methods

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include "data/DataWriterHDF5DataPointsTri3.hh" // USES DataWriterHDF5DataPointsTri3
#include "data/DataWriterHDF5DataPointsQuad4.hh" // USES DataWriterHDF5DataPointsQuad4
#include "data/DataWriterHDF5DataPointsTet4.hh" // USES DataWriterHDF5DataPointsTet4
#include "data/DataWriterHDF5DataPointsHex8.hh" // USES DataWriterHDF5DataPointsHex8


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5PointsTri3 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5PointsQuad4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5PointsTet4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5PointsHex8 );


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5PointsTri3::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5Points::setUp();
  _data = new DataWriterHDF5DataPointsTri3;
  _initialize();

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5PointsQuad4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5Points::setUp();
  _data = new DataWriterHDF5DataPointsQuad4;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5PointsTet4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5Points::setUp();
  _data = new DataWriterHDF5DataPointsTet4;
  _initialize();

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5PointsHex8::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterHDF5Points::setUp();
  _data = new DataWriterHDF5DataPointsHex8;
  _initialize();

  PYLITH_METHOD_END;
} // setUp

// End of file 
