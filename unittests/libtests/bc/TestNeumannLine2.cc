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

#include "TestNeumannLine2.hh" // Implementation of class methods

#include "data/NeumannDataLine2.hh" // USES NeumannDataLine2

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryPoint1D.hh" // USES GeometryPoint1D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestNeumannLine2 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestNeumannLine2::setUp(void)
{ // setUp
  TestNeumann::setUp();
  _data = new NeumannDataLine2();
  feassemble::GeometryPoint1D geometry;
  CPPUNIT_ASSERT(0 != _quadrature);
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
