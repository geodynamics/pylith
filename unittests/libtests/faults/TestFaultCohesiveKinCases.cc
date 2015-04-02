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

#include "TestFaultCohesiveKinCases.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature<Mesh>

#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D
#include "pylith/feassemble/GeometryTri3D.hh" // USES GeometryTri3D
#include "pylith/feassemble/GeometryQuad3D.hh" // USES GeometryQuad3D

// ----------------------------------------------------------------------
#include "data/CohesiveKinDataTri3.hh" // USES CohesiveKinDataTri3
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinTri3 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinTri3::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataTri3();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);
  
  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveKinDataTri3d.hh" // USES CohesiveKinDataTri3d
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinTri3d );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinTri3d::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataTri3d();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);
  
  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveKinDataTri3g.hh" // USES CohesiveKinDataTri3g
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinTri3g );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinTri3g::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataTri3g();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);
  
  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveKinDataQuad4.hh" // USES CohesiveKinDataQuad4
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinQuad4 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinQuad4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataQuad4();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveKinDataQuad4e.hh" // USES CohesiveKinDataQuad4e
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinQuad4e );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinQuad4e::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataQuad4e();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);
  
  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveKinDataQuad4i.hh" // USES CohesiveKinDataQuad4i
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinQuad4i );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinQuad4i::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataQuad4i();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);
  
  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveKinDataTet4.hh" // USES CohesiveKinDataTet4
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinTet4 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinTet4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataTet4();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryTri3D geometry;
  _quadrature->refGeometry(&geometry);

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveKinDataTet4e.hh" // USES CohesiveKinDataTet4e
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinTet4e );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinTet4e::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataTet4e();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryTri3D geometry;
  _quadrature->refGeometry(&geometry);

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveKinDataTet4f.hh" // USES CohesiveKinDataTet4f
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinTet4f );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinTet4f::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataTet4f();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryTri3D geometry;
  _quadrature->refGeometry(&geometry);

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveKinDataHex8.hh" // USES CohesiveKinDataHex8
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinHex8 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinHex8::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataHex8();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryQuad3D geometry;
  _quadrature->refGeometry(&geometry);
  
  PYLITH_METHOD_END;
} // setUp


// End of file 
