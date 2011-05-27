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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestExodusII.hh" // Implementation of class methods

#include "pylith/meshio/ExodusII.hh"

#include "pylith/utils/array.hh" // USES int_array, double_array, string_vector

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestExodusII );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestExodusII::testConstructor(void)
{ // testConstructor
  ExodusII exofile;
} // testConstructor

// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestExodusII::testFilename(void)
{ // testFilename
  ExodusII exofile;

  const char* filename = "hi.txt";
  exofile.filename(filename);
  CPPUNIT_ASSERT(0 == strcmp(filename, exofile.filename()));
} // testFilename

// ----------------------------------------------------------------------
// Test open() and close();
void
pylith::meshio::TestExodusII::testOpenClose(void)
{ // testOpenClose
  ExodusII exofile("data/twotri3.exo");
  CPPUNIT_ASSERT(exofile._file);

  exofile.close();
  CPPUNIT_ASSERT(!exofile._file);
} // testOpenClose

// ----------------------------------------------------------------------
// Test hasDim()
void
pylith::meshio::TestExodusII::testHasDim(void)
{ // testHasDim
  ExodusII exofile("data/twotri3.exo");

  CPPUNIT_ASSERT(exofile.hasDim("num_dim"));
  CPPUNIT_ASSERT(!exofile.hasDim("abcdefghijklm"));
} // testHasDim

// ----------------------------------------------------------------------
// Test hasAtt()
void
pylith::meshio::TestExodusII::testHasAtt(void)
{ // testHasAtt
  ExodusII exofile("data/twotri3.exo");

  CPPUNIT_ASSERT(exofile.hasAtt("version"));
  CPPUNIT_ASSERT(!exofile.hasAtt("abcdefghijklm"));
} // testHasAtt

// ----------------------------------------------------------------------
// Test hasVar()
void
pylith::meshio::TestExodusII::testHasVar(void)
{ // testHasVar
  ExodusII exofile("data/twotri3.exo");

  CPPUNIT_ASSERT(exofile.hasVar("connect1"));
  CPPUNIT_ASSERT(!exofile.hasVar("abcdefghijklm"));
} // testHasVar

// ----------------------------------------------------------------------
// Test getVar(double*).
void
pylith::meshio::TestExodusII::testGetVarDouble(void)
{ // testGetVarDouble
  const double coordsE[8] = { -1.0, 0.0, 0.0, 1.0,
			      0.0, -1.0, 1.0, 0.0 };

  const int ndims = 2;
  int dims[2];
  dims[0] = 2;
  dims[1] = 4;
  const int size = dims[0]*dims[1];
  double_array coords(size);

  ExodusII exofile("data/twotri3.exo");
  exofile.getVar(&coords[0], dims, ndims, "coord");

  const double tolerance = 1.0e-06;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(coordsE[i], coords[i], tolerance);
} // testGetVarDouble

// ----------------------------------------------------------------------
// Test getVar(int*).
void
pylith::meshio::TestExodusII::testGetVarInt(void)
{ // testGetVarDouble
  const int connectE[3] = { 3, 2, 4 };

  const int ndims = 2;
  int dims[2];
  dims[0] = 1;
  dims[1] = 3;
  const int size = dims[0]*dims[1];
  int_array connect(size);

  ExodusII exofile("data/twotri3.exo");
  exofile.getVar(&connect[0], dims, ndims, "connect2");

  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(connectE[i], connect[i]);
} // testGetVarDouble

// ----------------------------------------------------------------------
// Test getVar(string_vector).
void
pylith::meshio::TestExodusII::testGetVarString(void)
{ // testGetVarString
  const char* namesE[2] = { "x", "y" };

  const int dim = 2;
  string_vector names(2);

  ExodusII exofile("data/twotri3.exo");
  exofile.getVar(&names, dim, "coor_names");

  for (int i=0; i < dim; ++i)
    CPPUNIT_ASSERT_EQUAL(std::string(namesE[i]), names[i]);
} // testGetVarString


// End of file 
