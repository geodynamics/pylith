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

#include "TestFaultCohesiveKinQuad4e.hh" // Implementation of class methods

#include "data/CohesiveKinDataQuad4e.hh" // USES CohesiveKinDataQuad4e

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature<SubMesh>
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinQuad4e );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinQuad4e::setUp(void)
{ // setUp
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataQuad4e();

  CPPUNIT_ASSERT(0 != _quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);
  
  _flipFault = true;
} // setUp


// End of file 
