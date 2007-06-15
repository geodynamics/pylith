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

#include "TestFaultCohesiveKinTet4e.hh" // Implementation of class methods

#include "data/CohesiveKinDataTet4e.hh" // USES CohesiveKinDataTet4e

#include "pylith/feassemble/Quadrature2Din3D.hh" // USES Quadrature2Din3D
#include "pylith/feassemble/GeometryTri3D.hh" // USES GeometryTri3D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinTet4e );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinTet4e::setUp(void)
{ // setUp
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataTet4e();
  _quadrature = new feassemble::Quadrature2Din3D();
  CPPUNIT_ASSERT(0 != _quadrature);
  feassemble::GeometryTri3D geometry;
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
