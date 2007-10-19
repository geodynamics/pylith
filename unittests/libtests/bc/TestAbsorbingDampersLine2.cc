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

#include "TestAbsorbingDampersLine2.hh" // Implementation of class methods

#include "data/AbsorbingDampersDataLine2.hh" // USES AbsorbingDampersDataLine2

#include "pylith/feassemble/Quadrature0D.hh" // USES Quadrature1D
#include "pylith/feassemble/GeometryPoint1D.hh" // USES GeometryPoint1D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestAbsorbingDampersLine2 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestAbsorbingDampersLine2::setUp(void)
{ // setUp
  _data = new AbsorbingDampersDataLine2();
  _quadrature = new feassemble::Quadrature0D();
  CPPUNIT_ASSERT(0 != _quadrature);
  feassemble::GeometryPoint1D geometry;
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
