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

#include "TestFaultCohesiveDynCases.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature<Mesh>
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D
#include "pylith/feassemble/GeometryTri3D.hh" // USES GeometryTri3D
#include "pylith/feassemble/GeometryQuad3D.hh" // USES GeometryQuad3D

// ----------------------------------------------------------------------
#include "data/CohesiveDynDataTri3.hh" // USES CohesiveDynDataTri3
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveDynTri3 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveDynTri3::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestFaultCohesiveDyn::setUp();
  _data = new CohesiveDynDataTri3();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);
  
  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveDynDataTri3d.hh" // USES CohesiveDynDataTri3d
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveDynTri3d );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveDynTri3d::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestFaultCohesiveDyn::setUp();
  _data = new CohesiveDynDataTri3d();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);
  
  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveDynDataQuad4.hh" // USES CohesiveDynDataQuad4
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveDynQuad4 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveDynQuad4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestFaultCohesiveDyn::setUp();
  _data = new CohesiveDynDataQuad4();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);
  
  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveDynDataTet4.hh" // USES CohesiveDynDataTet4
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveDynTet4 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveDynTet4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestFaultCohesiveDyn::setUp();
  _data = new CohesiveDynDataTet4();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryTri3D geometry;
  _quadrature->refGeometry(&geometry);
  
  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveDynDataHex8.hh" // USES CohesiveDynLDataHex8
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveDynHex8 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveDynHex8::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestFaultCohesiveDyn::setUp();
  _data = new CohesiveDynDataHex8();

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryQuad3D geometry;
  _quadrature->refGeometry(&geometry);
  
  PYLITH_METHOD_END;
} // setUp


// End of file 
