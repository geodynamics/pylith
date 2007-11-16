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

#include "pylith/feassemble/Quadrature0D.hh" // USES Quadrature0D
#include "pylith/feassemble/GeometryPoint1D.hh" // USES GeometryPoint1D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestNeumannLine2 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestNeumannLine2::setUp(void)
{ // setUp
  _data = new NeumannDataLine2();
  _quadrature = new feassemble::Quadrature0D();
  CPPUNIT_ASSERT(0 != _quadrature);
  feassemble::GeometryPoint1D geometry;
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
