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

#include "TestFaultCohesiveKinTri3d.hh" // Implementation of class methods

#include "data/CohesiveKinDataTri3d.hh" // USES CohesiveKinDataTri3d

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature<SubMesh>
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinTri3d );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinTri3d::setUp(void)
{ // setUp
  TestFaultCohesiveKin::setUp();
  _data = new CohesiveKinDataTri3d();

  CPPUNIT_ASSERT(0 != _quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);
  
  _flipFault = true;
} // setUp


// End of file 
