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

#include "TestAbsorbingDampersTri3.hh" // Implementation of class methods

#include "data/AbsorbingDampersDataTri3.hh" // USES AbsorbingDampersDataTri3

#include "pylith/feassemble/Quadrature1Din2D.hh" // USES Quadrature1Din2D
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestAbsorbingDampersTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestAbsorbingDampersTri3::setUp(void)
{ // setUp
  _data = new AbsorbingDampersDataTri3();
  _quadrature = new feassemble::Quadrature1Din2D();
  CPPUNIT_ASSERT(0 != _quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
