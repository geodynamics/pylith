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

#include "TestFaultCohesiveKinSrcsLine2.hh" // Implementation of class methods

#include "data/CohesiveKinSrcsDataLine2.hh" // USES CohesiveKinDataSrcsLine2

#include "pylith/feassemble/Quadrature0D.hh" // USES Quadrature0D
#include "pylith/feassemble/GeometryPoint1D.hh" // USES GeometryPoint1D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinSrcsLine2 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinSrcsLine2::setUp(void)
{ // setUp
  TestFaultCohesiveKinSrcs::setUp();
  _data = new CohesiveKinSrcsDataLine2();
  _quadrature = new feassemble::Quadrature0D();
  CPPUNIT_ASSERT(0 != _quadrature);
  feassemble::GeometryPoint1D geometry;
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
