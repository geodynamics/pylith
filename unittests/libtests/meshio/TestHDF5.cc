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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestHDF5.hh" // Implementation of class methods

#include "pylith/meshio/HDF5.hh" // USES HDF5

#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 8
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

  HDF5 two("test.h5", H5F_ACC_TRUNC);
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

  h5.open("test.h5", H5F_ACC_TRUNC);
  CPPUNIT_ASSERT(h5._file >= 0);
  CPPUNIT_ASSERT(h5.isOpen());
  
  h5.close();
  CPPUNIT_ASSERT(-1 == h5._file);
  CPPUNIT_ASSERT(!h5.isOpen());

  h5.open("test.h5", H5F_ACC_RDONLY);
  CPPUNIT_ASSERT(h5._file >= 0);
  CPPUNIT_ASSERT(h5.isOpen());
  h5.close();
  CPPUNIT_ASSERT(-1 == h5._file);
  CPPUNIT_ASSERT(!h5.isOpen());
} // testOpenClose

// ----------------------------------------------------------------------
// Test hasGroup()
void
pylith::meshio::TestHDF5::testHasGroup(void)
{ // testHasGroup
  HDF5 h5("test.h5", H5F_ACC_TRUNC);

  h5.createGroup("/mygroup");
  h5.close();

  h5.open("test.h5", H5F_ACC_RDONLY);
  CPPUNIT_ASSERT(h5.hasGroup("/mygroup"));
  CPPUNIT_ASSERT(!h5.hasGroup("/notmygroup"));
  h5.close();
} // testHasGroup

// ----------------------------------------------------------------------
// Test hasDataset()
void
pylith::meshio::TestHDF5::testHasDataset(void)
{ // testHasDataset
  HDF5 h5("test.h5", H5F_ACC_TRUNC);

  const hsize_t ndims = 1;
  const hsize_t dims[ndims] = { 2 };
  h5.createDataset("/", "data", dims, dims, ndims, H5T_NATIVE_INT);
  h5.close();

  h5.open("test.h5", H5F_ACC_RDONLY);
  CPPUNIT_ASSERT(h5.hasDataset("/data"));
  CPPUNIT_ASSERT(!h5.hasDataset("/nodata"));
  h5.close();
} // testHasDataset

// ----------------------------------------------------------------------
// Test getDatasetDims().
void
pylith::meshio::TestHDF5::testGetDatasetDims(void)
{ // testGetDatasetDims
  HDF5 h5("test.h5", H5F_ACC_TRUNC);

  const int ndimsE = 2;
  const hsize_t dimsE[ndimsE] = { 3, 2 };
  const hsize_t dimsChunkE[ndimsE] = { 1, 2 };
  h5.createDataset("/", "data", dimsE, dimsChunkE, ndimsE, H5T_NATIVE_INT);
  h5.close();

  h5.open("test.h5", H5F_ACC_RDONLY);
  hsize_t* dims = 0;
  int ndims = 0;
  h5.getDatasetDims(&dims, &ndims, "/", "data");
  h5.close();
  CPPUNIT_ASSERT_EQUAL(ndimsE, ndims);

  for (int i=0; i < ndimsE; ++i)
    CPPUNIT_ASSERT_EQUAL(dimsE[i], dims[i]);

  delete[] dims; dims = 0;
} // testGetDatasetDims

// ----------------------------------------------------------------------
// Test getGroupDatasets().
void
pylith::meshio::TestHDF5::testGetGroupDatasets(void)
{ // testGetGroupDatasets
  const int ngroupsE = 3;
  const char* namesE[3] = { "dataA",
			    "dataB",
			    "dataC" };
  const hsize_t ndims = 2;
  const hsize_t dims[ndims] = { 3, 2 };
  const hsize_t dimsChunk[ndims] = { 1, 2 };

  HDF5 h5("test.h5", H5F_ACC_TRUNC);
  h5.createGroup("/mygroup");
  for (int i=0; i < ngroupsE; ++i)
    h5.createDataset("/mygroup", namesE[i], dims, dimsChunk, ndims,
		     H5T_NATIVE_INT);
  h5.close();

  string_vector names;
  h5.open("test.h5", H5F_ACC_RDONLY);
  h5.getGroupDatasets(&names, "/mygroup");
  h5.close();

  const int ngroups = names.size();
  CPPUNIT_ASSERT_EQUAL(ngroupsE, ngroups);
  for (int i=0; i < ngroups; ++i)
    CPPUNIT_ASSERT_EQUAL(std::string(namesE[i]), names[i]);
} // testGetGroupDatasets

// ----------------------------------------------------------------------
// Test createGroup()
void
pylith::meshio::TestHDF5::testCreateGroup(void)
{ // testCreateGroup
  HDF5 h5("test.h5", H5F_ACC_TRUNC);

  h5.createGroup("/mygroup");
  h5.close();

  h5.open("test.h5", H5F_ACC_RDONLY);
#if defined(PYLITH_HDF5_USE_API_18)
  hid_t group = H5Gopen2(h5._file, "/mygroup", H5P_DEFAULT);
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
  HDF5 h5("test.h5", H5F_ACC_TRUNC);

  const hsize_t ndims = 1;
  const hsize_t dims[ndims] = { 2 };
  h5.createDataset("/", "data", dims, dims, ndims, H5T_NATIVE_INT);

  const PylithScalar scalarE = 2.5;
  h5.writeAttribute("/data", "myscalar", (void*)&scalarE, H5T_NATIVE_DOUBLE);
  h5.close();

  const PylithScalar tolerance = 1.0e-06;
  h5.open("test.h5", H5F_ACC_RDONLY);
  PylithScalar scalar = 0;
  h5.readAttribute("/data", "myscalar", (void*)&scalar, H5T_NATIVE_DOUBLE);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(scalarE, scalar, tolerance);
  h5.close();
} // testAttributeScalar

// ----------------------------------------------------------------------
// Test writeAttribute(string) and readAttribute(string).
void
pylith::meshio::TestHDF5::testAttributeString(void)
{ // testAttributeString
  HDF5 h5("test.h5", H5F_ACC_TRUNC);

  const hsize_t ndims = 1;
  const hsize_t dims[ndims] = { 2 };
  h5.createDataset("/", "data", dims, dims, ndims, H5T_NATIVE_INT);

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
  HDF5 h5("test.h5", H5F_ACC_TRUNC);

  const hsize_t ndims = 2;
  const hsize_t dims[ndims] = { 3, 2 };
  const hsize_t dimsChunk[ndims] = { 1, 2 };
  h5.createDataset("/", "data", dims, dimsChunk, ndims, H5T_NATIVE_INT);
  h5.close();

  h5.open("test.h5", H5F_ACC_RDONLY);
#if defined(PYLITH_HDF5_USE_API_18)
  hid_t group = H5Gopen2(h5._file, "/", H5P_DEFAULT);
#else
  hid_t group = H5Gopen(h5._file, "/");
#endif
  CPPUNIT_ASSERT(group >= 0);
#if defined(PYLITH_HDF5_USE_API_18)
  hid_t dataset = H5Dopen2(group, "data", H5P_DEFAULT);
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
// Test writeDatasetChunk() and readDatasetChunk().
void
pylith::meshio::TestHDF5::testDatasetChunk(void)
{ // testDatasetChunk
  const int ndimsE = 3;
  const hsize_t dimsE[ndimsE] = { 4, 2, 3 };
  const hsize_t dimsChunkE[ndimsE] = { 1, 2, 3 };

  // Create data.
  hsize_t nitems = dimsE[0];
  hsize_t nitemsS = 1;
  for (int i=1; i < ndimsE; ++i) {
    nitems *= dimsE[i];
    nitemsS *= dimsE[i];
  } // for
  const hsize_t sizeBytes = nitems * H5Tget_size(H5T_NATIVE_INT);
  int* valuesE = (nitems > 0) ? new int[nitems] : 0;
  for (int i=0; i < nitems; ++i)
    valuesE[i] = 2 * i + 1;

  HDF5 h5("test.h5", H5F_ACC_TRUNC);
  h5.createDataset("/", "data", dimsE, dimsChunkE, ndimsE, H5T_NATIVE_INT);

  for (int i=0; i < dimsE[0]; ++i)
    h5.writeDatasetChunk("/", "data", (void*)&valuesE[i*nitemsS], 
			 dimsE, dimsChunkE, ndimsE, i, H5T_NATIVE_INT);
  h5.close();
  
  int ndims = 0;
  hsize_t* dims = 0;
  int* values = 0;
  h5.open("test.h5", H5F_ACC_RDONLY);
  for (int i=0; i < dimsE[0]; ++i) {
    h5.readDatasetChunk("/", "data", (char**)&values, &dims, &ndims, i, 
			H5T_NATIVE_INT);
    CPPUNIT_ASSERT_EQUAL(ndimsE, ndims);
    for (int iDim=1; iDim < ndims; ++iDim)
      CPPUNIT_ASSERT_EQUAL(dimsE[iDim], dims[iDim]);

    for (int ii=0; ii < nitemsS; ++ii)
      CPPUNIT_ASSERT_EQUAL(valuesE[i*nitemsS+ii], values[ii]);
  } // for

  delete[] values; values = 0;
  delete[] dims; dims = 0;
  delete[] valuesE; valuesE = 0;

  h5.close();
} // testDatasetChunk

// ----------------------------------------------------------------------
// Test createDatasetRawExternal() and updateDatasetRawExternal().
void
pylith::meshio::TestHDF5::testDatasetRawExternal(void)
{ // testDatasetRawExternal
  const int ndimsE = 2;
  const hsize_t dimsE[ndimsE] = { 6, 3 };
  hsize_t dims[ndimsE];

  // Create raw data file
  hsize_t nitems = dimsE[0];
  for (int i=1; i < ndimsE; ++i)
    nitems *= dimsE[i];
  const hsize_t sizeBytes = nitems * H5Tget_size(H5T_NATIVE_INT);
  int* valuesE = (nitems > 0) ? new int[nitems] : 0;
  for (int i=0; i < nitems; ++i)
    valuesE[i] = 2 * i + 1;
  std::ofstream fout("test.dat");
  fout.write((char*)valuesE, sizeBytes);
  fout.close();

  HDF5 h5("test.h5", H5F_ACC_TRUNC);
  dims[0] = H5S_UNLIMITED;
  dims[1] = dimsE[1];
  h5.createDatasetRawExternal("/", "data", "test.dat", dims, ndimsE,
			      H5T_NATIVE_INT);
  h5.extendDatasetRawExternal("/", "data", dimsE, ndimsE);
  h5.close();

  h5.open("test.h5", H5F_ACC_RDONLY);
#if defined(PYLITH_HDF5_USE_API_18)
  hid_t group = H5Gopen2(h5._file, "/", H5P_DEFAULT);
#else
  hid_t group = H5Gopen(h5._file, "/");
#endif
  CPPUNIT_ASSERT(group >= 0);
#if defined(PYLITH_HDF5_USE_API_18)
  hid_t dataset = H5Dopen2(group, "data", H5P_DEFAULT);
#else
  hid_t dataset = H5Dopen(group, "data");
#endif
  CPPUNIT_ASSERT(dataset >= 0);

  hid_t dataspace = H5Dget_space(dataset);
  CPPUNIT_ASSERT(dataspace >= 0);
  
  const int ndims = H5Sget_simple_extent_ndims(dataspace);
  CPPUNIT_ASSERT_EQUAL(ndimsE, ndims);
  H5Sget_simple_extent_dims(dataspace, dims, 0);
  for (int i=0; i < ndims; ++i)
    CPPUNIT_ASSERT_EQUAL(dimsE[i], dims[i]);

  int* values = (nitems > 0) ? new int[nitems] : 0;
  herr_t err = H5Dread(dataset, H5T_NATIVE_INT, dataspace, dataspace, 
		       H5P_DEFAULT, (void*)values);
  CPPUNIT_ASSERT(err >= 0);

  for (int i=0; i < nitems; ++i)
    CPPUNIT_ASSERT_EQUAL(valuesE[i], values[i]);
  delete[] valuesE; valuesE = 0;
  delete[] values; values = 0;

  err = H5Sclose(dataspace);
  CPPUNIT_ASSERT(err >= 0);
  err = H5Dclose(dataset);
  CPPUNIT_ASSERT(err >= 0);
  err = H5Gclose(group);
  CPPUNIT_ASSERT(err >= 0);
  h5.close();
} // testDatasetRawExternal


// End of file 
