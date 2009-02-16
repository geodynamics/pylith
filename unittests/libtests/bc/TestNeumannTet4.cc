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

#include "TestNeumannTet4.hh" // Implementation of class methods

#include "data/NeumannDataTet4.hh" // USES NeumannDataTet4

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryTri3D.hh" // USES GeometryTri3D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestNeumannTet4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestNeumannTet4::setUp(void)
{ // setUp
  TestNeumann::setUp();
  _data = new NeumannDataTet4();
  feassemble::GeometryTri3D geometry;
  CPPUNIT_ASSERT(0 != _quadrature);
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
