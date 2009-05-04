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

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryQuad3D.hh" // USES GeometryQuad3D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestNeumannHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestNeumannHex8::setUp(void)
{ // setUp
  TestNeumann::setUp();
  _data = new NeumannDataHex8();
  feassemble::GeometryQuad3D geometry;
  CPPUNIT_ASSERT(0 != _quadrature);
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
