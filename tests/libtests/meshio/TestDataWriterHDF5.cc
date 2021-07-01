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

#include "TestDataWriterHDF5.hh" // Implementation of class methods

#include "pylith/utils/types.hh" // HASA PylithScalar
#include "pylith/utils/error.h" // HASA PYLITH_METHOD_BEGIN/END

#include <hdf5.h> // USES HDF5 API

#if H5_VERSION_GE(1,12,0)
#define PYLITH_HDF5_USE_API_112
#endif

// ------------------------------------------------------------------------------------------------
herr_t
pylith_meshio_TestDataWriterHDF5_checkObject(hid_t id,
                                             const char* name,
                                             const H5O_info_t* info,
                                             void* data) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(info);
    CPPUNIT_ASSERT(data);

    hid_t* file = (hid_t*) data;CPPUNIT_ASSERT(H5Iis_valid(*file));
    herr_t err = 0;

    switch (info->type) {
    case H5O_TYPE_GROUP: {
        hid_t group = H5Gopen2(*file, name, H5P_DEFAULT);CPPUNIT_ASSERT(group >= 0);
        err = H5Gclose(group);CPPUNIT_ASSERT(err >= 0);
        break;
    } // group
    case H5O_TYPE_DATASET: {
        // Get expected dataset.
        hid_t datasetE = H5Dopen2(id, name, H5P_DEFAULT);CPPUNIT_ASSERT(datasetE >= 0);
        hid_t dataspaceE = H5Dget_space(datasetE);CPPUNIT_ASSERT(dataspaceE >= 0);
        const int ndimsE = H5Sget_simple_extent_ndims(dataspaceE);CPPUNIT_ASSERT(ndimsE > 0);
        hsize_t* dimsE = (ndimsE > 0) ? new hsize_t[ndimsE] : 0;
        const int ndimsECheck = H5Sget_simple_extent_dims(dataspaceE, dimsE, 0);CPPUNIT_ASSERT_EQUAL(ndimsE, ndimsECheck);
        int sizeE = (ndimsE > 0 && dimsE[0] > 0) ? 1 : 0;
        for (int i = 0; i < ndimsE; ++i) {
            sizeE *= dimsE[i];
        } // for

        // Get dataset
        hid_t dataset = H5Dopen2(*file, name, H5P_DEFAULT);CPPUNIT_ASSERT_MESSAGE(name, dataset >= 0);
        hid_t dataspace = H5Dget_space(dataset);CPPUNIT_ASSERT(dataspace >= 0);
        const int ndims = H5Sget_simple_extent_ndims(dataspace);CPPUNIT_ASSERT(ndims > 0);
        hsize_t* dims = (ndims > 0) ? new hsize_t[ndims] : 0;
        const int ndimsCheck = H5Sget_simple_extent_dims(dataspace, dims, 0);CPPUNIT_ASSERT_EQUAL(ndims, ndimsCheck);
        int size = (ndims > 0 && dims[0] > 0) ? 1 : 0;
        for (int i = 0; i < ndims; ++i) {
            size *= dims[i];
        } // for

        // Check dimensions.
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in number of dimensions.", ndimsE, ndims);
        for (int i = 0; i < ndimsE; ++i) {
            std::ostringstream msg;
            msg << "Mismatch in dimension "<< i << " of dataset '" << name << "'." << std::endl;
            CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str().c_str(), dimsE[i], dims[i]);
        } // for

        // Check the expected datatype
        hid_t datatypeE = H5Dget_type(datasetE);CPPUNIT_ASSERT(datatypeE >= 0);
        hid_t dataclassE = H5Tget_class(datatypeE);CPPUNIT_ASSERT(dataclassE >= 0);

        hid_t datatype = H5Dget_type(dataset);CPPUNIT_ASSERT(datatype >= 0);
        hid_t dataclass = H5Tget_class(datatype);CPPUNIT_ASSERT(dataclass >= 0);
        CPPUNIT_ASSERT_EQUAL(dataclassE, dataclass);

        switch (dataclassE) {
        case H5T_FLOAT: {
            double* dataE = (sizeE > 0) ? new double[sizeE] : 0;CPPUNIT_ASSERT(sizeE > 0);
            err = H5Dread(datasetE, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*) dataE);CPPUNIT_ASSERT(err >= 0);

            double* data = (size > 0) ? new double[size] : 0;CPPUNIT_ASSERT(size > 0);
            err = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*) data);CPPUNIT_ASSERT(err >= 0);

            CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in dataset total size.", sizeE, size);

            // Compare data values.
            const double tolerance = 1.0e-6;
            for (int i = 0; i < size; ++i) {
                const double toleranceV = std::max(tolerance, tolerance*dataE[i]);
                std::ostringstream msg;
                msg << "Mismatch in dataset '" << name << "' at index " << i << ".";
                CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(msg.str().c_str(), dataE[i], double(data[i]), toleranceV);
            } // for

            delete[] dataE;dataE = 0;
            delete[] data;data = 0;

            break;
        } // H5T_DOUBLE

        case H5T_STRING: {
            const int slenE = H5Tget_size(datatypeE);CPPUNIT_ASSERT(slenE > 0);
            sizeE *= slenE;

            const int slen = H5Tget_size(datatype);CPPUNIT_ASSERT(slen > 0);
            size *= slen;

            CPPUNIT_ASSERT_EQUAL(slenE, slen);
            CPPUNIT_ASSERT_EQUAL(sizeE, size);

            char* dataE = (sizeE > 0) ? new char[sizeE] : 0;CPPUNIT_ASSERT(sizeE > 0);
            err = H5Dread(datasetE, datatypeE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataE);CPPUNIT_ASSERT(err >= 0);

            char* data = (size > 0) ? new char[size] : 0;CPPUNIT_ASSERT(size > 0);
            err = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);CPPUNIT_ASSERT(err >= 0);

            for (int i = 0; i < size; ++i) {
                std::ostringstream msg;
                msg << "Mismatch in dataset '" << name << "' at index " << i << ".";
                CPPUNIT_ASSERT_EQUAL(dataE[i], data[i]);
            } // for

            delete[] dataE;dataE = 0;
            delete[] data;data = 0;

            break;
        } // H5T_C_S1

        default:
            CPPUNIT_ASSERT(false);
        } // switch

        err = H5Sclose(dataspaceE);CPPUNIT_ASSERT(err >= 0);
        err = H5Dclose(datasetE);CPPUNIT_ASSERT(err >= 0);

        err = H5Sclose(dataspace);CPPUNIT_ASSERT(err >= 0);
        err = H5Dclose(dataset);CPPUNIT_ASSERT(err >= 0);

        delete[] dimsE;dimsE = 0;
        delete[] dims;dims = 0;

        break;
    } // dataset
    default:
        CPPUNIT_ASSERT(false);
    } // switch

    PYLITH_METHOD_RETURN(0);
} // checkObject


// ------------------------------------------------------------------------------------------------
// Check HDF5 file against archived file.
void
pylith::meshio::TestDataWriterHDF5::checkFile(const char* filename) {
    PYLITH_METHOD_BEGIN;

    const std::string filenameE = "data/" + std::string(filename);

    herr_t err = 0;

    hid_t fileE = H5Fopen(filenameE.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);CPPUNIT_ASSERT(fileE >= 0);

    hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);CPPUNIT_ASSERT(file >= 0);
#if defined(PYLITH_HDF5_USE_API_112)
    // Traverse recursively file with expected values.
    err = H5Ovisit(fileE, H5_INDEX_NAME, H5_ITER_NATIVE, pylith_meshio_TestDataWriterHDF5_checkObject, (void*) &file, H5O_INFO_ALL);CPPUNIT_ASSERT(err >= 0);
#else
    err = H5Ovisit(fileE, H5_INDEX_NAME, H5_ITER_NATIVE, pylith_meshio_TestDataWriterHDF5_checkObject, (void*) &file);CPPUNIT_ASSERT(err >= 0);
#endif
    err = H5Fclose(fileE);
    CPPUNIT_ASSERT(err >= 0);

    err = H5Fclose(file);
    CPPUNIT_ASSERT(err >= 0);

    PYLITH_METHOD_END;
} // checkFile


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::TestDataWriterHDF5_Data::TestDataWriterHDF5_Data(void) :
    opencloseFilename(NULL),
    vertexFilename(NULL),
    cellFilename(NULL) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::TestDataWriterHDF5_Data::~TestDataWriterHDF5_Data(void) {}


// End of file
