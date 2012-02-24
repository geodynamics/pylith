// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestAbsorbingDampersLine2.hh" // Implementation of class methods

#include "data/AbsorbingDampersDataLine2.hh" // USES AbsorbingDampersDataLine2

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryPoint1D.hh" // USES GeometryPoint0D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestAbsorbingDampersLine2 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestAbsorbingDampersLine2::setUp(void)
{ // setUp
  TestAbsorbingDampers::setUp();
  _data = new AbsorbingDampersDataLine2();
  feassemble::GeometryPoint1D geometry;
  CPPUNIT_ASSERT(0 != _quadrature);
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
