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

#include "TestFaultCohesiveDynLTet4.hh" // Implementation of class methods

#include "data/CohesiveDynLDataTet4.hh" // USES CohesiveDynLDataTet4

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature<SubMesh>
#include "pylith/feassemble/GeometryTri3D.hh" // USES GeometryTri3D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveDynLTet4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveDynLTet4::setUp(void)
{ // setUp
  TestFaultCohesiveDynL::setUp();
  _data = new CohesiveDynLDataTet4();

  CPPUNIT_ASSERT(0 != _quadrature);
  feassemble::GeometryTri3D geometry;
  _quadrature->refGeometry(&geometry);
  
  _flipFault = true;
} // setUp


// End of file 
