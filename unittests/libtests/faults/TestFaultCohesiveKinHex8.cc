// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFaultCohesiveKinHex8.hh" // Implementation of class methods

#include "data/CohesiveKinDataHex8.hh" // USES CohesiveKinDataHex8

#include "pylith/feassemble/Quadrature2Din3D.hh" // USES Quadrature2Din3D
#include "pylith/feassemble/GeometryQuad3D.hh" // USES GeometryQuad3D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinHex8::setUp(void)
{ // setUp
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataHex8();
  _quadrature = new feassemble::Quadrature2Din3D();
  CPPUNIT_ASSERT(0 != _quadrature);
  feassemble::GeometryQuad3D geometry;
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
