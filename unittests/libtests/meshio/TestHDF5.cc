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

#include "TestHDF5.hh" // Implementation of class methods

#include "pylith/meshio/HDF5.hh" // USES HDF5

#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 8 && !defined(H5_USE_API_16)
#define PYLITH_HDF5_USE_API_18
#endif

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestHDF5 );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::meshio::TestHDF5::testConstructor(void)
{ // testConstructor
  HDF5 one;
  CPPUNIT_ASSERT(-1 == one._file);

  HDF5 two("test.h5", H5F_ACC_TRUNC, true);
  CPPUNIT_ASSERT(two._file >= 0);
  two.close();

  HDF5 three("test.h5", H5F_ACC_RDONLY);
  CPPUNIT_ASSERT(three._file >= 0);
} // testConstructor

// ----------------------------------------------------------------------
// Test open() and close().
void
pylith::meshio::TestHDF5::testOpenClose(void)
{ // testOpenClose
  HDF5 h5;
  CPPUNIT_ASSERT(-1 == h5._file);

  h5.open("test.h5", H5F_ACC_TRUNC, true);
  CPPUNIT_ASSERT(h5._file >= 0);
  
  h5.close();
  CPPUNIT_ASSERT(-1 == h5._file);

  h5.open("test.h5", H5F_ACC_RDONLY);
  CPPUNIT_ASSERT(h5._file >= 0);
  h5.close();
  CPPUNIT_ASSERT(-1 == h5._file);
} // testOpenClose

// ----------------------------------------------------------------------
// Test createGroup()
void
pylith::meshio::TestHDF5::testCreateGroup(void)
{ // testCreateGroup
  HDF5 h5("test.h5", H5F_ACC_TRUNC, true);

  h5.createGroup("/mygroup");
  h5.close();

  h5.open("test.h5", H5F_ACC_RDONLY);
#if defined(PYLITH_HDF5_USE_API_18)
  hid_t group = H5Gopen(h5._file, "/mygroup", H5P_DEFAULT);
#else
  hid_t group = H5Gopen(h5._file, "/mygroup");
#endif
  CPPUNIT_ASSERT(group >= 0);
  herr_t err = H5Gclose(group);
  CPPUNIT_ASSERT(err >= 0);
  h5.close();
} // testCreateGroup

// ----------------------------------------------------------------------
// Test writeAttribute(scalar) and readAttribute(scalar).
void
pylith::meshio::TestHDF5::testAttributeScalar(void)
{ // testAttributeScalar
  HDF5 h5("test.h5", H5F_ACC_TRUNC, true);

  const hsize_t ndims = 1;
  const hsize_t dims[ndims] = { 2 };
  h5.createDataset("/", "data", dims, ndims, H5T_NATIVE_INT);

  const double scalarE = 2.5;
  h5.writeAttribute("/data", "myscalar", (void*)&scalarE, H5T_NATIVE_DOUBLE);
  h5.close();

  const double tolerance = 1.0e-06;
  h5.open("test.h5", H5F_ACC_RDONLY);
  double scalar = 0;
  h5.readAttribute("/data", "myscalar", (void*)&scalar, H5T_NATIVE_DOUBLE);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(scalarE, scalar, tolerance);
  h5.close();
} // testAttributeScalar

// ----------------------------------------------------------------------
// Test writeAttribute(string) and readAttribute(string).
void
pylith::meshio::TestHDF5::testAttributeString(void)
{ // testAttributeString
  HDF5 h5("test.h5", H5F_ACC_TRUNC, true);

  const hsize_t ndims = 1;
  const hsize_t dims[ndims] = { 2 };
  h5.createDataset("/", "data", dims, ndims, H5T_NATIVE_INT);

  const std::string valueE = "abcd";
  h5.writeAttribute("/data", "mystring", valueE.c_str());
  h5.close();

  h5.open("test.h5", H5F_ACC_RDONLY);
  std::string value = h5.readAttribute("/data", "mystring");
  CPPUNIT_ASSERT_EQUAL(valueE, value);
  h5.close();
} // testAttributeString

// ----------------------------------------------------------------------
// Test createDataset().
void
pylith::meshio::TestHDF5::testCreateDataset(void)
{ // testCreateDataset
  HDF5 h5("test.h5", H5F_ACC_TRUNC, true);

  const hsize_t ndims = 1;
  const hsize_t dims[ndims] = { 2 };
  h5.createDataset("/", "data", dims, ndims, H5T_NATIVE_INT);
  h5.close();

  h5.open("test.h5", H5F_ACC_RDONLY);
#if defined(PYLITH_HDF5_USE_API_18)
  hid_t group = H5Gopen(h5._file, "/", H5P_DEFAULT);
#else
  hid_t group = H5Gopen(h5._file, "/");
#endif
  CPPUNIT_ASSERT(group >= 0);
#if defined(PYLITH_HDF5_USE_API_18)
  hid_t dataset = H5Dopen(group, "data", H5P_DEFAULT);
#else
  hid_t dataset = H5Dopen(group, "data");
#endif
  CPPUNIT_ASSERT(dataset >= 0);
  herr_t err = H5Dclose(dataset);
  CPPUNIT_ASSERT(err >= 0);
  err = H5Gclose(group);
  CPPUNIT_ASSERT(err >= 0);
  h5.close();
} // testCreateDataset

// ----------------------------------------------------------------------
// Test createDatasetRawExternal.
void
pylith::meshio::TestHDF5::testCreateDatasetRawExternal(void)
{ // testCreateDatasetRawExternal
  const hsize_t ndims = 1;
  const hsize_t dims[ndims] = { 6 };

  // Create raw data file
  hsize_t nitems = 0;
  for (int i=0; i < ndims; ++i)
    nitems += dims[i];
  const hsize_t sizeBytes = nitems * H5Tget_size(H5T_NATIVE_INT);
  int* values = (nitems > 0) ? new int[nitems] : 0;
  for (int i=0; i < nitems; ++i)
    values[i] = 2 * i + 1;
  std::ofstream fout("test.dat");
  fout.write((char*)values, sizeBytes);
  fout.close();
  delete[] values; values = 0;

  HDF5 h5("test.h5", H5F_ACC_TRUNC, true);
  h5.createDatasetRawExternal("/", "data", "test.dat", dims, ndims, H5T_NATIVE_INT);
  h5.close();

  h5.open("test.h5", H5F_ACC_RDONLY);
#if defined(PYLITH_HDF5_USE_API_18)
  hid_t group = H5Gopen(h5._file, "/", H5P_DEFAULT);
#else
  hid_t group = H5Gopen(h5._file, "/");
#endif
  CPPUNIT_ASSERT(group >= 0);
#if defined(PYLITH_HDF5_USE_API_18)
  hid_t dataset = H5Dopen(group, "data", H5P_DEFAULT);
#else
  hid_t dataset = H5Dopen(group, "data");
#endif
  CPPUNIT_ASSERT(dataset >= 0);
  herr_t err = H5Dclose(dataset);
  CPPUNIT_ASSERT(err >= 0);
  err = H5Gclose(group);
  CPPUNIT_ASSERT(err >= 0);
  h5.close();
} // testCreateDatasetRawExternal

// ----------------------------------------------------------------------
// Test writeDatasetSlice.
void
pylith::meshio::TestHDF5::testWriteDatasetSlice(void)
{ // testWriteDatasetSlice
} // testWriteDatasetSlice


// End of file 
