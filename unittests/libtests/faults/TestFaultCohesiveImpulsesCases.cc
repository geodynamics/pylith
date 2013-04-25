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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFaultCohesiveImpulsesCases.hh" // Implementation of class methods

#include "data/CohesiveImpulsesDataTri3.hh"
#include "data/CohesiveImpulsesDataQuad4.hh"
#include "data/CohesiveImpulsesDataTet4.hh"
#include "data/CohesiveImpulsesDataHex8.hh"

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature<SubMesh>
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D
#include "pylith/feassemble/GeometryTri3D.hh" // USES GeometryTri3D
#include "pylith/feassemble/GeometryQuad3D.hh" // USES GeometryQuad3D


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveImpulsesTri3 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveImpulsesQuad4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveImpulsesTet4 );
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveImpulsesHex8 );


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveImpulsesTri3::setUp(void)
{ // setUp
  TestFaultCohesiveImpulses::setUp();
  _data = new CohesiveImpulsesDataTri3;

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry); 

  _flipFault = true;
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveImpulsesQuad4::setUp(void)
{ // setUp
  TestFaultCohesiveImpulses::setUp();
  _data = new CohesiveImpulsesDataQuad4;

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry); 
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveImpulsesTet4::setUp(void)
{ // setUp
  TestFaultCohesiveImpulses::setUp();
  _data = new CohesiveImpulsesDataTet4;

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryTri3D geometry;
  _quadrature->refGeometry(&geometry); 
} // setUp


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveImpulsesHex8::setUp(void)
{ // setUp
  TestFaultCohesiveImpulses::setUp();
  _data = new CohesiveImpulsesDataHex8;

  CPPUNIT_ASSERT(_quadrature);
  feassemble::GeometryQuad3D geometry;
  _quadrature->refGeometry(&geometry); 

  _flipFault = true;
} // setUp


// End of file 
