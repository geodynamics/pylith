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

#include "TestNeumannQuad4.hh" // Implementation of class methods

#include "data/NeumannDataQuad4.hh" // USES NeumannDataQuad4

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestNeumannQuad4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestNeumannQuad4::setUp(void)
{ // setUp
  TestNeumann::setUp();
  _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(0 != _quadrature);
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
