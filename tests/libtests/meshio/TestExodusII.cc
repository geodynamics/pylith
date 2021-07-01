// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>
#include <stdexcept>

#include "TestExodusII.hh" // Implementation of class methods

#include "pylith/meshio/ExodusII.hh"

#include "pylith/utils/array.hh" // USES int_array, scalar_array, string_vector
#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestExodusII );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestExodusII::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  ExodusII exofile;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestExodusII::testFilename(void)
{ // testFilename
  PYLITH_METHOD_BEGIN;

  ExodusII exofile;

  const char* filename = "hi.txt";
  exofile.filename(filename);
  CPPUNIT_ASSERT(0 == strcmp(filename, exofile.filename()));

  PYLITH_METHOD_END;
} // testFilename

// ----------------------------------------------------------------------
// Test open() and close();
void
pylith::meshio::TestExodusII::testOpenClose(void)
{ // testOpenClose
  PYLITH_METHOD_BEGIN;

  ExodusII exofile("data/twotri3_13.0.exo");
  CPPUNIT_ASSERT(exofile._file);

  exofile.close();
  CPPUNIT_ASSERT(!exofile._file);

  // Attempt to open file that doesn't exist.
  exofile.filename("fail.exo");
  CPPUNIT_ASSERT_THROW(exofile.open(), std::runtime_error);

  // Attempt to close file with bad handle.
  exofile._file = 1;
  CPPUNIT_ASSERT_THROW(exofile.close(), std::runtime_error);

  PYLITH_METHOD_END;
} // testOpenClose

// ----------------------------------------------------------------------
// Test hasDim()
void
pylith::meshio::TestExodusII::testHasDim(void)
{ // testHasDim
  PYLITH_METHOD_BEGIN;

  ExodusII exofile("data/twotri3_13.0.exo");

  int id = -1;
  CPPUNIT_ASSERT(exofile.hasDim("num_dim", &id));
  CPPUNIT_ASSERT(id >= 0);
  CPPUNIT_ASSERT(!exofile.hasDim("abcdefghijklm", &id));
  CPPUNIT_ASSERT(-1 == id);

  PYLITH_METHOD_END;
} // testHasDim

// ----------------------------------------------------------------------
// Test hasAtt()
void
pylith::meshio::TestExodusII::testHasAtt(void)
{ // testHasAtt
  PYLITH_METHOD_BEGIN;

  ExodusII exofile("data/twotri3_13.0.exo");

  int id = -1;
  CPPUNIT_ASSERT(exofile.hasAtt("version", &id));
  CPPUNIT_ASSERT(id >= 0);
  CPPUNIT_ASSERT(!exofile.hasAtt("abcdefghijklm", &id));
  CPPUNIT_ASSERT(-1 == id);

  PYLITH_METHOD_END;
} // testHasAtt

// ----------------------------------------------------------------------
// Test hasVar()
void
pylith::meshio::TestExodusII::testHasVar(void)
{ // testHasVar
  PYLITH_METHOD_BEGIN;

  ExodusII exofile("data/twotri3_13.0.exo");

  int id = -1;
  CPPUNIT_ASSERT(exofile.hasVar("connect1", &id));
  CPPUNIT_ASSERT(id >= 0);
  CPPUNIT_ASSERT(!exofile.hasVar("abcdefghijklm", &id));
  CPPUNIT_ASSERT(-1 == id);

  PYLITH_METHOD_END;
} // testHasVar

// ----------------------------------------------------------------------
// Test getVar(PylithScalar*).
void
pylith::meshio::TestExodusII::testGetVarDouble(void)
{ // testGetVarDouble
  PYLITH_METHOD_BEGIN;

  const PylithScalar coordsE[8] = { -1.0, 0.0, 0.0, 1.0,
			      0.0, -1.0, 1.0, 0.0 };

  const int ndims = 2;
  int dims[2];
  dims[0] = 2;
  dims[1] = 4;
  const int size = dims[0]*dims[1];
  scalar_array coords(size);

  ExodusII exofile("data/twotri3_12.2.exo");
  exofile.getVar(&coords[0], dims, ndims, "coord");

  const PylithScalar tolerance = 1.0e-06;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(coordsE[i], coords[i], tolerance);

  // Attempt to get variable that doesn't exist.
  CPPUNIT_ASSERT_THROW(exofile.getVar(&coords[0], dims, ndims, "aabbcc"), std::runtime_error);

  // Attempt to get variable with wrong number of dimensions.
  CPPUNIT_ASSERT_THROW(exofile.getVar(&coords[0], dims, ndims+1, "coord"), std::runtime_error);

  // Attempt to get variable with wrong dimension.
  dims[0] = 99;
  CPPUNIT_ASSERT_THROW(exofile.getVar(&coords[0], dims, ndims, "coord"), std::runtime_error);

  PYLITH_METHOD_END;
} // testGetVarDouble

// ----------------------------------------------------------------------
// Test getVar(int*).
void
pylith::meshio::TestExodusII::testGetVarInt(void)
{ // testGetVarDouble
  PYLITH_METHOD_BEGIN;

  const int connectE[3] = { 3, 2, 4 };

  const int ndims = 2;
  int dims[2];
  dims[0] = 1;
  dims[1] = 3;
  const int size = dims[0]*dims[1];
  int_array connect(size);

  ExodusII exofile("data/twotri3_13.0.exo");
  exofile.getVar(&connect[0], dims, ndims, "connect2");

  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(connectE[i], connect[i]);

  // Attempt to get variable that doesn't exist.
  CPPUNIT_ASSERT_THROW(exofile.getVar(&connect[0], dims, ndims, "aabbcc"), std::runtime_error);

  // Attempt to get variable with wrong number of dimensions.
  CPPUNIT_ASSERT_THROW(exofile.getVar(&connect[0], dims, ndims+1, "connect2"), std::runtime_error);

  // Attempt to get variable with wrong dimension.
  dims[0] = 99;
  CPPUNIT_ASSERT_THROW(exofile.getVar(&connect[0], dims, ndims, "connect2"), std::runtime_error);

  PYLITH_METHOD_END;
} // testGetVarDouble

// ----------------------------------------------------------------------
// Test getVar(string_vector).
void
pylith::meshio::TestExodusII::testGetVarString(void)
{ // testGetVarString
  PYLITH_METHOD_BEGIN;

  const char* namesE[2] = { "x", "y" };

  const int dim = 2;
  string_vector names(2);

  ExodusII exofile("data/twotri3_13.0.exo");
  exofile.getVar(&names, dim, "coor_names");

  for (int i=0; i < dim; ++i)
    CPPUNIT_ASSERT_EQUAL(std::string(namesE[i]), names[i]);

  // Attempt to get variable that doesn't exist.
  CPPUNIT_ASSERT_THROW(exofile.getVar(&names, dim, "aabbcc"), std::runtime_error);

  // Attempt to get variable with wrong number of dimensions.
  CPPUNIT_ASSERT_THROW(exofile.getVar(&names, dim+1, "coord_names"), std::runtime_error);

  PYLITH_METHOD_END;
} // testGetVarString


// End of file
