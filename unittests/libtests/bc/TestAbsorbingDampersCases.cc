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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestAbsorbingDampersCases.hh" // Implementation of cases

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestAbsorbingDampersTri3 );

#include "data/AbsorbingDampersDataTri3.hh" // USES AbsorbingDampersDataTri3
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

// Setup testing data.
void
pylith::bc::TestAbsorbingDampersTri3::setUp(void)
{ // setUp
  TestAbsorbingDampers::setUp();
  _data = new AbsorbingDampersDataTri3();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestAbsorbingDampersQuad4 );

#include "data/AbsorbingDampersDataQuad4.hh" // USES AbsorbingDampersDataQuad4
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

// Setup testing data.
void
pylith::bc::TestAbsorbingDampersQuad4::setUp(void)
{ // setUp
  TestAbsorbingDampers::setUp();
  _data = new AbsorbingDampersDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestAbsorbingDampersTet4 );

#include "data/AbsorbingDampersDataTet4.hh" // USES AbsorbingDampersDataTet4
#include "pylith/feassemble/GeometryTri3D.hh" // USES GeometryTri3D

// Setup testing data.
void
pylith::bc::TestAbsorbingDampersTet4::setUp(void)
{ // setUp
  TestAbsorbingDampers::setUp();
  _data = new AbsorbingDampersDataTet4();
  feassemble::GeometryTri3D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestAbsorbingDampersHex8 );

#include "data/AbsorbingDampersDataHex8.hh" // USES AbsorbingDampersDataHex8
#include "pylith/feassemble/GeometryQuad3D.hh" // USES GeometryQuad3D

// Setup testing data.
void
pylith::bc::TestAbsorbingDampersHex8::setUp(void)
{ // setUp
  TestAbsorbingDampers::setUp();
  _data = new AbsorbingDampersDataHex8();
  feassemble::GeometryQuad3D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
