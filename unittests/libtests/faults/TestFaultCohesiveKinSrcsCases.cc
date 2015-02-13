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

#include "TestFaultCohesiveKinSrcsCases.hh" // Implementation of class methods


#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature<Mesh>

#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D
#include "pylith/feassemble/GeometryTri3D.hh" // USES GeometryTri3D
#include "pylith/feassemble/GeometryQuad3D.hh" // USES GeometryQuad3D

// ----------------------------------------------------------------------
#include "data/CohesiveKinSrcsDataTri3.hh" // USES CohesiveKinSrcsDataTri3
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinSrcsTri3 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinSrcsTri3::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestFaultCohesiveKinSrcs::setUp();
  _data = new CohesiveKinSrcsDataTri3();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);
  
  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveKinSrcsDataQuad4.hh" // USES CohesiveKinDataQuad4
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinSrcsQuad4 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinSrcsQuad4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestFaultCohesiveKinSrcs::setUp();
  _data = new CohesiveKinSrcsDataQuad4();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveKinSrcsDataTet4.hh" // USES CohesiveKinDataTet4
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinSrcsTet4 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinSrcsTet4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestFaultCohesiveKinSrcs::setUp();
  _data = new CohesiveKinSrcsDataTet4();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryTri3D geometry;
  _quadrature->refGeometry(&geometry);

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveKinSrcsDataHex8.hh" // USES CohesiveKinDataHex8
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinSrcsHex8 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinSrcsHex8::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestFaultCohesiveKinSrcs::setUp();
  _data = new CohesiveKinSrcsDataHex8();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryQuad3D geometry;
  _quadrature->refGeometry(&geometry);
  
  PYLITH_METHOD_END;
} // setUp


// End of file 
