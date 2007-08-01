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

#include "TestNeumannHex8.hh" // Implementation of class methods

#include "data/NeumannDataHex8.hh" // USES NeumannDataHex8

#include "pylith/feassemble/Quadrature2Din3D.hh" // USES Quadrature2Din3D
#include "pylith/feassemble/GeometryQuad3D.hh" // USES GeometryQuad3D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestNeumannHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestNeumannHex8::setUp(void)
{ // setUp
  _data = new NeumannDataHex8();
  _quadrature = new feassemble::Quadrature2Din3D();
  CPPUNIT_ASSERT(0 != _quadrature);
  feassemble::GeometryQuad3D geometry;
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
