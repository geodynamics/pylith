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

#include "TestNeumannCases.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestNeumannTri3 );

#include "data/NeumannDataTri3.hh" // USES NeumannDataTri3
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

// Setup testing data.
void
pylith::bc::TestNeumannTri3::setUp(void)
{ // setUp
  TestNeumann::setUp();
  _data = new NeumannDataTri3();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestNeumannQuad4 );

#include "data/NeumannDataQuad4.hh" // USES NeumannDataQuad4
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

// Setup testing data.
void
pylith::bc::TestNeumannQuad4::setUp(void)
{ // setUp
  TestNeumann::setUp();
  _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestNeumannTet4 );

#include "data/NeumannDataTet4.hh" // USES NeumannDataTet4
#include "pylith/feassemble/GeometryTri3D.hh" // USES GeometryTri3D

// Setup testing data.
void
pylith::bc::TestNeumannTet4::setUp(void)
{ // setUp
  TestNeumann::setUp();
  _data = new NeumannDataTet4();
  feassemble::GeometryTri3D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestNeumannHex8 );

#include "data/NeumannDataHex8.hh" // USES NeumannDataHex8
#include "pylith/feassemble/GeometryQuad3D.hh" // USES GeometryQuad3D

// Setup testing data.
void
pylith::bc::TestNeumannHex8::setUp(void)
{ // setUp
  TestNeumann::setUp();
  _data = new NeumannDataHex8();
  feassemble::GeometryQuad3D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);
} // setUp


// End of file 
