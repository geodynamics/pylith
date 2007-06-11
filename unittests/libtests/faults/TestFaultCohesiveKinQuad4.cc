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

#include "TestFaultCohesiveKinQuad4.hh" // Implementation of class methods

#include "data/CohesiveKinDataQuad4.hh" // USES CohesiveKinDataQuad4

#include "pylith/feassemble/Quadrature1Din2D.hh" // USES Quadrature1Din2D
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinQuad4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinQuad4::setUp(void)
{ // setUp
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataQuad4();
  _quadrature = new feassemble::Quadrature1Din2D();
  CPPUNIT_ASSERT(0 != _quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
