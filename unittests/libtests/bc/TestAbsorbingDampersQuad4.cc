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

#include "TestAbsorbingDampersQuad4.hh" // Implementation of class methods

#include "data/AbsorbingDampersDataQuad4.hh" // USES AbsorbingDampersDataQuad4

#include "pylith/feassemble/Quadrature1Din2D.hh" // USES Quadrature1Din2D
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestAbsorbingDampersQuad4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestAbsorbingDampersQuad4::setUp(void)
{ // setUp
  _data = new AbsorbingDampersDataQuad4();
  _quadrature = new feassemble::Quadrature1Din2D();
  CPPUNIT_ASSERT(0 != _quadrature);
  feassemble::GeometryLine2D geometry;
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
