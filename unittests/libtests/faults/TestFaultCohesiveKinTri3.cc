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

#include "TestFaultCohesiveKinTri3.hh" // Implementation of class methods

#include "data/CohesiveKinDataTri3.hh" // USES CohesiveKinDataTri3

#include "pylith/feassemble/Quadrature1D.hh" // USES Quadrature1D
#include "pylith/feassemble/GeometryLine1D.hh" // USES GeometryLine1D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinTri3::setUp(void)
{ // setUp
  _data = new CohesiveKinDataTri3();
  _quadrature = new feassemble::Quadrature1D();
  CPPUNIT_ASSERT(0 != _quadrature);
  feassemble::GeometryLine1D geometry;
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
