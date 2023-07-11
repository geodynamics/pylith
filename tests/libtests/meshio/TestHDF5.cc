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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/meshio/HDF5.hh" // USES HDF5
#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <fstream> // USES ofstream

#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 8
#define PYLITH_HDF5_USE_API_18
#endif

namespace pylith {
    namespace meshio {
        class TestHDF5;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestHDF5 : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Test constructor.
    static
    void testConstructor(void);

    /// Test open() and close().
    static
    void testOpenClose(void);

    /// Test hasGroup().
    static
    void testHasGroup(void);

    /// Test hasDataset().
    static
    void testHasDataset(void);

    /// Test getDatasetDims().
    static
    void testGetDatasetDims(void);

    /// Test getGroupDatasets().
    static
    void testGetGroupDatasets(void);

    /// Test createGroup()
    static
    void testCreateGroup(void);

    /// Test writeAttribute(scalar) and readAttribute(scalar).
    static
    void testAttributeScalar(void);

    /// Test createDataset().
    static
    void testCreateDataset(void);

    /// Test writeDatasetChunk() and readDatasetChunk().
    static
    void testDatasetChunk(void);

    /// Test createDatasetRawExternal() and updateDatasetRawExternal().
    static
    void testDatasetRawExternal(void);

    /// Test writeAttribute(string) and readAttribute(string).
    static
    void testAttributeString(void);

    /// Test writeDataset(string) and readDataset(string).
    static
    void testDatasetString(void);

}; // class TestHDF5

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestHDF5::testConstructor", "[TestHDF5]") {
    pylith::meshio::TestHDF5::testConstructor();
}
TEST_CASE("TestHDF5::testOpenClose", "[TestHDF5]") {
    pylith::meshio::TestHDF5::testOpenClose();
}
TEST_CASE("TestHDF5::testHasGroup", "[TestHDF5]") {
    pylith::meshio::TestHDF5::testHasGroup();
}
TEST_CASE("TestHDF5::testHasDataset", "[TestHDF5]") {
    pylith::meshio::TestHDF5::testHasDataset();
}
TEST_CASE("TestHDF5::testGetDatasetDims", "[TestHDF5]") {
    pylith::meshio::TestHDF5::testGetDatasetDims();
}
TEST_CASE("TestHDF5::testGetGroupDatasets", "[TestHDF5]") {
    pylith::meshio::TestHDF5::testGetGroupDatasets();
}
TEST_CASE("TestHDF5::testCreateGroup", "[TestHDF5]") {
    pylith::meshio::TestHDF5::testCreateGroup();
}
TEST_CASE("TestHDF5::testAttributeScalar", "[TestHDF5]") {
    pylith::meshio::TestHDF5::testAttributeScalar();
}
TEST_CASE("TestHDF5::testCreateDataset", "[TestHDF5]") {
    pylith::meshio::TestHDF5::testCreateDataset();
}
TEST_CASE("TestHDF5::testDatasetChunk", "[TestHDF5]") {
    pylith::meshio::TestHDF5::testDatasetChunk();
}
TEST_CASE("TestHDF5::testDatasetRawExternal", "[TestHDF5]") {
    pylith::meshio::TestHDF5::testDatasetRawExternal();
}
TEST_CASE("TestHDF5::testAttributeString", "[TestHDF5]") {
    pylith::meshio::TestHDF5::testAttributeString();
}
TEST_CASE("TestHDF5::testDatasetString", "[TestHDF5]") {
    pylith::meshio::TestHDF5::testDatasetString();
}

// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::meshio::TestHDF5::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    HDF5 one;
    CHECK(-1 == one._file);

    HDF5 two("test.h5", H5F_ACC_TRUNC);
    CHECK(two._file >= 0);
    two.close();

    HDF5 three("test.h5", H5F_ACC_RDONLY);
    CHECK(three._file >= 0);

    PYLITH_METHOD_END;
} // testConstructor


// ------------------------------------------------------------------------------------------------
// Test open() and close().
void
pylith::meshio::TestHDF5::testOpenClose(void) {
    PYLITH_METHOD_BEGIN;

    HDF5 h5;
    CHECK(-1 == h5._file);

    h5.open("test.h5", H5F_ACC_TRUNC);
    CHECK(h5._file >= 0);
    CHECK(h5.isOpen());

    h5.close();
    CHECK(-1 == h5._file);
    CHECK(!h5.isOpen());

    h5.open("test.h5", H5F_ACC_RDONLY);
    CHECK(h5._file >= 0);
    CHECK(h5.isOpen());
    h5.close();
    CHECK(-1 == h5._file);
    CHECK(!h5.isOpen());

    PYLITH_METHOD_END;
} // testOpenClose


// ------------------------------------------------------------------------------------------------
// Test hasGroup()
void
pylith::meshio::TestHDF5::testHasGroup(void) {
    PYLITH_METHOD_BEGIN;

    HDF5 h5("test.h5", H5F_ACC_TRUNC);

    h5.createGroup("/mygroup");
    h5.close();

    h5.open("test.h5", H5F_ACC_RDONLY);
    CHECK(h5.hasGroup("/mygroup"));
    CHECK(!h5.hasGroup("/notmygroup"));
    h5.close();

    PYLITH_METHOD_END;
} // testHasGroup


// ------------------------------------------------------------------------------------------------
// Test hasDataset()
void
pylith::meshio::TestHDF5::testHasDataset(void) {
    PYLITH_METHOD_BEGIN;

    HDF5 h5("test.h5", H5F_ACC_TRUNC);

    const hsize_t ndims = 1;
    const hsize_t dims[ndims] = { 2 };
    h5.createDataset("/", "data", dims, dims, ndims, H5T_NATIVE_INT);
    h5.close();

    h5.open("test.h5", H5F_ACC_RDONLY);
    CHECK(h5.hasDataset("/data"));
    CHECK(!h5.hasDataset("/nodata"));
    h5.close();

    PYLITH_METHOD_END;
} // testHasDataset


// ------------------------------------------------------------------------------------------------
// Test getDatasetDims().
void
pylith::meshio::TestHDF5::testGetDatasetDims(void) {
    PYLITH_METHOD_BEGIN;

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
    REQUIRE(ndimsE == ndims);

    for (int i = 0; i < ndimsE; ++i) {
        CHECK(dimsE[i] == dims[i]);
    }

    delete[] dims;dims = 0;

    PYLITH_METHOD_END;
} // testGetDatasetDims


// ------------------------------------------------------------------------------------------------
// Test getGroupDatasets().
void
pylith::meshio::TestHDF5::testGetGroupDatasets(void) {
    PYLITH_METHOD_BEGIN;

    const int ngroupsE = 3;
    const char* namesE[3] = {
        "dataA",
        "dataB",
        "dataC",
    };
    const hsize_t ndims = 2;
    const hsize_t dims[ndims] = { 3, 2 };
    const hsize_t dimsChunk[ndims] = { 1, 2 };

    HDF5 h5("test.h5", H5F_ACC_TRUNC);
    h5.createGroup("/mygroup");
    for (int i = 0; i < ngroupsE; ++i) {
        h5.createDataset("/mygroup", namesE[i], dims, dimsChunk, ndims,
                         H5T_NATIVE_INT);
    }
    h5.close();

    string_vector names;
    h5.open("test.h5", H5F_ACC_RDONLY);
    h5.getGroupDatasets(&names, "/mygroup");
    h5.close();

    const int ngroups = names.size();
    REQUIRE(ngroupsE == ngroups);
    for (int i = 0; i < ngroups; ++i) {
        CHECK(std::string(namesE[i]) == names[i]);
    }

    PYLITH_METHOD_END;
} // testGetGroupDatasets


// ------------------------------------------------------------------------------------------------
// Test createGroup()
void
pylith::meshio::TestHDF5::testCreateGroup(void) {
    PYLITH_METHOD_BEGIN;

    HDF5 h5("test.h5", H5F_ACC_TRUNC);

    h5.createGroup("/mygroup");
    h5.close();

    h5.open("test.h5", H5F_ACC_RDONLY);
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t group = H5Gopen2(h5._file, "/mygroup", H5P_DEFAULT);
#else
    hid_t group = H5Gopen(h5._file, "/mygroup");
#endif
    CHECK(group >= 0);
    herr_t err = H5Gclose(group);
    CHECK(err >= 0);
    h5.close();

    PYLITH_METHOD_END;
} // testCreateGroup


// ------------------------------------------------------------------------------------------------
// Test writeAttribute(scalar) and readAttribute(scalar).
void
pylith::meshio::TestHDF5::testAttributeScalar(void) {
    PYLITH_METHOD_BEGIN;

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
    CHECK_THAT(scalar, Catch::Matchers::WithinAbs(scalarE, tolerance));
    h5.close();

    PYLITH_METHOD_END;
} // testAttributeScalar


// ------------------------------------------------------------------------------------------------
// Test createDataset().
void
pylith::meshio::TestHDF5::testCreateDataset(void) {
    PYLITH_METHOD_BEGIN;

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
    CHECK(group >= 0);
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen2(group, "data", H5P_DEFAULT);
#else
    hid_t dataset = H5Dopen(group, "data");
#endif
    CHECK(dataset >= 0);
    herr_t err = H5Dclose(dataset);
    CHECK(err >= 0);
    err = H5Gclose(group);
    CHECK(err >= 0);
    h5.close();

    PYLITH_METHOD_END;
} // testCreateDataset


// ------------------------------------------------------------------------------------------------
// Test writeDatasetChunk() and readDatasetChunk().
void
pylith::meshio::TestHDF5::testDatasetChunk(void) {
    PYLITH_METHOD_BEGIN;

    const int ndimsE = 3;
    const hsize_t dimsE[ndimsE] = { 4, 2, 3 };
    const hsize_t dimsChunkE[ndimsE] = { 1, 2, 3 };

    // Create data.
    hsize_t nitems = dimsE[0];
    hsize_t nitemsS = 1;
    for (int i = 1; i < ndimsE; ++i) {
        nitems *= dimsE[i];
        nitemsS *= dimsE[i];
    } // for
    int* valuesE = (nitems > 0) ? new int[nitems] : 0;
    for (size_t i = 0; i < nitems; ++i) {
        valuesE[i] = 2 * i + 1;
    }

    HDF5 h5("test.h5", H5F_ACC_TRUNC);
    h5.createDataset("/", "data", dimsE, dimsChunkE, ndimsE, H5T_NATIVE_INT);

    for (size_t i = 0; i < dimsE[0]; ++i) {
        h5.writeDatasetChunk("/", "data", (void*)&valuesE[i*nitemsS],
                             dimsE, dimsChunkE, ndimsE, i, H5T_NATIVE_INT);
    }
    h5.close();

    int ndims = 0;
    hsize_t* dims = 0;
    int* values = 0;
    h5.open("test.h5", H5F_ACC_RDONLY);
    for (size_t i = 0; i < dimsE[0]; ++i) {
        h5.readDatasetChunk("/", "data", (char**)&values, &dims, &ndims, i,
                            H5T_NATIVE_INT);
        REQUIRE(ndimsE == ndims);
        for (int iDim = 1; iDim < ndims; ++iDim) {
            CHECK(dimsE[iDim] == dims[iDim]);
        }

        for (size_t ii = 0; ii < nitemsS; ++ii) {
            CHECK(valuesE[i*nitemsS+ii] == values[ii]);
        }
    } // for

    delete[] values;values = 0;
    delete[] dims;dims = 0;
    delete[] valuesE;valuesE = 0;

    h5.close();

    PYLITH_METHOD_END;
} // testDatasetChunk


// ------------------------------------------------------------------------------------------------
// Test createDatasetRawExternal() and updateDatasetRawExternal().
void
pylith::meshio::TestHDF5::testDatasetRawExternal(void) {
    PYLITH_METHOD_BEGIN;

    const int ndimsE = 2;
    const hsize_t dimsE[ndimsE] = { 6, 3 };
    hsize_t dims[ndimsE];

    // Create raw data file
    hsize_t nitems = dimsE[0];
    for (int i = 1; i < ndimsE; ++i) {
        nitems *= dimsE[i];
    }
    const hsize_t sizeBytes = nitems * H5Tget_size(H5T_NATIVE_INT);
    int* valuesE = (nitems > 0) ? new int[nitems] : 0;
    for (size_t i = 0; i < nitems; ++i) {
        valuesE[i] = 2 * i + 1;
    }
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
    CHECK(group >= 0);
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen2(group, "data", H5P_DEFAULT);
#else
    hid_t dataset = H5Dopen(group, "data");
#endif
    CHECK(dataset >= 0);

    hid_t dataspace = H5Dget_space(dataset);
    CHECK(dataspace >= 0);

    const int ndims = H5Sget_simple_extent_ndims(dataspace);
    REQUIRE(ndimsE == ndims);
    H5Sget_simple_extent_dims(dataspace, dims, 0);
    for (int i = 0; i < ndims; ++i) {
        CHECK(dimsE[i] == dims[i]);
    }

    int* values = (nitems > 0) ? new int[nitems] : 0;
    herr_t err = H5Dread(dataset, H5T_NATIVE_INT, dataspace, dataspace,
                         H5P_DEFAULT, (void*)values);
    CHECK(err >= 0);

    for (size_t i = 0; i < nitems; ++i) {
        CHECK(valuesE[i] == values[i]);
    }
    delete[] valuesE;valuesE = 0;
    delete[] values;values = 0;

    err = H5Sclose(dataspace);
    CHECK(err >= 0);
    err = H5Dclose(dataset);
    CHECK(err >= 0);
    err = H5Gclose(group);
    CHECK(err >= 0);
    h5.close();

    PYLITH_METHOD_END;
} // testDatasetRawExternal


// ------------------------------------------------------------------------------------------------
// Test writeAttribute(string) and readAttribute(string).
void
pylith::meshio::TestHDF5::testAttributeString(void) {
    PYLITH_METHOD_BEGIN;

    HDF5 h5("test.h5", H5F_ACC_TRUNC);

    const hsize_t ndims = 1;
    const hsize_t dims[ndims] = { 2 };
    h5.createDataset("/", "data", dims, dims, ndims, H5T_NATIVE_INT);

    const std::string valueE = "abcd";
    h5.writeAttribute("/data", "mystring", valueE.c_str());
    h5.close();

    h5.open("test.h5", H5F_ACC_RDONLY);
    std::string value = h5.readAttribute("/data", "mystring");
    CHECK(valueE == value);
    h5.close();

    PYLITH_METHOD_END;
} // testAttributeString


// ------------------------------------------------------------------------------------------------
// Test writeDataset(string) and readDataset(string).
void
pylith::meshio::TestHDF5::testDatasetString(void) {
    PYLITH_METHOD_BEGIN;

    HDF5 h5("test.h5", H5F_ACC_TRUNC);

    const size_t nstrings = 3;
    const char* dataE[nstrings] = {"abc", "defg", "hijkl" };

    h5.writeDataset("/", "data", dataE, nstrings);
    h5.close();

    h5.open("test.h5", H5F_ACC_RDONLY);
    const pylith::string_vector& data = h5.readDataset("/", "data");
    h5.close();

    REQUIRE(nstrings == data.size());
    for (size_t i = 0; i < nstrings; ++i) {
        const std::string& stringE = dataE[i];
        CHECK(stringE == data[i]);
    } // for

    PYLITH_METHOD_END;
} // testDatasetString


// End of file
