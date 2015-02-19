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

#include "TestFaultCohesiveImpulsesCases.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature<Mesh>

#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D
#include "pylith/feassemble/GeometryTri3D.hh" // USES GeometryTri3D
#include "pylith/feassemble/GeometryQuad3D.hh" // USES GeometryQuad3D


// ----------------------------------------------------------------------
#include "data/CohesiveImpulsesDataTri3.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveImpulsesTri3 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveImpulsesTri3::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestFaultCohesiveImpulses::setUp();
  _data = new CohesiveImpulsesDataTri3;

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry); 

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveImpulsesDataQuad4.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveImpulsesQuad4 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveImpulsesQuad4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestFaultCohesiveImpulses::setUp();
  _data = new CohesiveImpulsesDataQuad4;

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry); 

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveImpulsesDataTet4.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveImpulsesTet4 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveImpulsesTet4::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestFaultCohesiveImpulses::setUp();
  _data = new CohesiveImpulsesDataTet4;

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryTri3D geometry;
  _quadrature->refGeometry(&geometry); 

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/CohesiveImpulsesDataHex8.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveImpulsesHex8 );

// Setup testing data.
void
pylith::faults::TestFaultCohesiveImpulsesHex8::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestFaultCohesiveImpulses::setUp();
  _data = new CohesiveImpulsesDataHex8;

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryQuad3D geometry;
  _quadrature->refGeometry(&geometry); 

  PYLITH_METHOD_END;
} // setUp


// End of file 
