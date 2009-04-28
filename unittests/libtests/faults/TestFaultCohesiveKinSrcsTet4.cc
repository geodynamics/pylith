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

#include "TestFaultCohesiveKinSrcsTet4.hh" // Implementation of class methods

#include "data/CohesiveKinSrcsDataTet4.hh" // USES CohesiveKinDataTet4

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature<SubMesh>
#include "pylith/feassemble/GeometryTri3D.hh" // USES GeometryTri3D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinSrcsTet4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinSrcsTet4::setUp(void)
{ // setUp
  TestFaultCohesiveKinSrcs::setUp();
  _data = new CohesiveKinSrcsDataTet4();

  CPPUNIT_ASSERT(0 != _quadrature);
  feassemble::GeometryTri3D geometry;
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
