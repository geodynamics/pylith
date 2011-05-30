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

#include "TestDataWriterHDF5.hh" // Implementation of class methods

#include <hdf5.h> // USES HDF5 API

#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 8
#define PYLITH_HDF5_USE_API_18
#endif

// ----------------------------------------------------------------------
herr_t
pylith_meshio_TestDataWriterHDF5_checkObject(hid_t id, 
					     const char* name, 
					     const H5O_info_t* info,
					     void* data)
{ // checkObject
  CPPUNIT_ASSERT(info);  
  CPPUNIT_ASSERT(data);


  hid_t* file = (hid_t*) data;
  CPPUNIT_ASSERT(H5Iis_valid(*file));
  herr_t err = 0;

  switch (info->type) {
  case H5O_TYPE_GROUP : {
    hid_t group = H5Gopen2(*file, name, H5P_DEFAULT);
    CPPUNIT_ASSERT(group >= 0);
    err = H5Gclose(group);
    CPPUNIT_ASSERT(err >= 0);
    break;
  } // group
  case H5O_TYPE_DATASET : {
    // Get expected dataset.
    hid_t datasetE = H5Dopen2(id, name, H5P_DEFAULT);
    CPPUNIT_ASSERT(datasetE >= 0);
    hid_t dataspaceE = H5Dget_space(datasetE);
    CPPUNIT_ASSERT(dataspaceE >= 0);
    const int ndimsE = H5Sget_simple_extent_ndims(dataspaceE);
    CPPUNIT_ASSERT(ndimsE > 0);
    hsize_t* dimsE = (ndimsE > 0) ? new hsize_t[ndimsE] : 0;
    const int ndimsECheck = H5Sget_simple_extent_dims(dataspaceE, dimsE, 0);
    CPPUNIT_ASSERT_EQUAL(ndimsE, ndimsECheck);
    int sizeE = (ndimsE > 0 && dimsE[0] > 0) ? 1 : 0;
    for (int i=0; i < ndimsE; ++i)
      sizeE *= dimsE[i];
    double* dataE = (sizeE > 0) ? new double[sizeE] : 0;
    CPPUNIT_ASSERT(sizeE > 0);
    err = H5Dread(datasetE, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		  H5P_DEFAULT, (void*) dataE);
    CPPUNIT_ASSERT(err >= 0);

    err = H5Sclose(dataspaceE);
    CPPUNIT_ASSERT(err >= 0);
    err = H5Dclose(datasetE);
    CPPUNIT_ASSERT(err >= 0);
    
    // Get test dataset.
    hid_t dataset = H5Dopen2(*file, name, H5P_DEFAULT);
    CPPUNIT_ASSERT(dataset >= 0);
    hid_t dataspace = H5Dget_space(dataset);
    CPPUNIT_ASSERT(dataspace >= 0);

    const int ndims = H5Sget_simple_extent_ndims(dataspace);
    CPPUNIT_ASSERT(ndims > 0);
    hsize_t* dims = (ndims > 0) ? new hsize_t[ndims] : 0;
    const int ndimsCheck = H5Sget_simple_extent_dims(dataspace, dims, 0);
    CPPUNIT_ASSERT_EQUAL(ndims, ndimsCheck);
    int size = (ndims > 0 && dims[0] > 0) ? 1 : 0;
    for (int i=0; i < ndims; ++i)
      size *= dims[i];
    double* data = (size > 0) ? new double[size] : 0;
    CPPUNIT_ASSERT(size > 0);
    err = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		  H5P_DEFAULT, (void*) data);
    CPPUNIT_ASSERT(err >= 0);

    err = H5Sclose(dataspace);
    CPPUNIT_ASSERT(err >= 0);
    err = H5Dclose(dataset);
    CPPUNIT_ASSERT(err >= 0);

    // Compare dimensions.
    CPPUNIT_ASSERT_EQUAL(ndimsE, ndims);
    for (int i=0; i < ndimsE; ++i)
      CPPUNIT_ASSERT_EQUAL(dimsE[i], dims[i]);

    // Compare data values.
    const double tolerance = 1.0e-6;
    CPPUNIT_ASSERT_EQUAL(sizeE, size);
    for (int i=0; i < size; ++i)
      if (dataE[i] != 0.0)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, data[i]/dataE[i], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(dataE[i], data[i], tolerance);

    break;
  } // dataset
  default :
    CPPUNIT_ASSERT(false);
  } // switch

  return 0;
} // checkObject

// ----------------------------------------------------------------------
// Check HDF5 file against archived file.
void
pylith::meshio::TestDataWriterHDF5::checkFile(const char* filename)
{ // checkFile

  const std::string filenameE = "data/" + std::string(filename);

  herr_t err = 0;

  hid_t fileE = H5Fopen(filenameE.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  CPPUNIT_ASSERT(fileE >= 0);

  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  CPPUNIT_ASSERT(file >= 0);
#if defined(PYLITH_HDF5_USE_API_18)

  // Traverse recursively file with expected values.
  err = H5Ovisit(fileE, H5_INDEX_NAME, H5_ITER_NATIVE, 
		 pylith_meshio_TestDataWriterHDF5_checkObject, (void*) &file);
  CPPUNIT_ASSERT(err >= 0);

#else
#endif
  err = H5Fclose(fileE);
  CPPUNIT_ASSERT(err >= 0);

  err = H5Fclose(file);
  CPPUNIT_ASSERT(err >= 0);

} // checkFile


// End of file 
