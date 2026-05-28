// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/meshio/HDF5.hh" // USES HDF5
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <memory> // USGS std::unique_ptr
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

    h5.createDataset<int>("/", "data", HDF5::DatasetShape{{2}, {2}});
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

    const size_t ndimsE = 2;
    const hsize_t dimsE[ndimsE] = { 3, 2 };
    h5.createDataset<int>("/", "data", HDF5::DatasetShape{{3, 2}, {1, 2}});
    h5.close();

    h5.open("test.h5", H5F_ACC_RDONLY);
    HDF5::DatasetDims result = h5.getDatasetDims("/", "data");
    h5.close();
    REQUIRE(ndimsE == result.dims.size());

    for (size_t i = 0; i < ndimsE; ++i) {
        CHECK(dimsE[i] == result.dims[i]);
    } // for

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

    HDF5 h5("test.h5", H5F_ACC_TRUNC);
    h5.createGroup("/mygroup");
    for (int i = 0; i < ngroupsE; ++i) {
        h5.createDataset<int>("/mygroup", namesE[i], HDF5::DatasetShape{{3, 2}, {1, 2}});
    }
    h5.close();

    h5.open("test.h5", H5F_ACC_RDONLY);
    pylith::string_vector names = h5.getGroupDatasets("/mygroup");
    h5.close();

    const int ngroups = names.size();
    REQUIRE(ngroupsE == ngroups);
    for (int i = 0; i < ngroups; ++i) {
        CHECK(std::string(namesE[i]) == names[i]);
    } // for

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

    _HDF5::H5Handle group(-1, H5Gclose);

    h5.open("test.h5", H5F_ACC_RDONLY);
#if defined(PYLITH_HDF5_USE_API_18)
    group.id = H5Gopen2(h5._file, "/mygroup", H5P_DEFAULT);
#else
    group.id = H5Gopen(h5._file, "/mygroup");
#endif
    h5.close();
    CHECK(group.id >= 0);

    PYLITH_METHOD_END;
} // testCreateGroup


// ------------------------------------------------------------------------------------------------
// Test writeAttribute(scalar) and readAttribute(scalar).
void
pylith::meshio::TestHDF5::testAttributeScalar(void) {
    PYLITH_METHOD_BEGIN;

    HDF5 h5("test.h5", H5F_ACC_TRUNC);

    h5.createDataset<int>("/", "data", HDF5::DatasetShape{{2}, {2}});

    const PylithScalar scalarE = 2.5;
    h5.writeAttribute<PylithScalar>("/data", "myscalar", scalarE);
    h5.close();

    const PylithScalar tolerance = 1.0e-06;
    h5.open("test.h5", H5F_ACC_RDONLY);
    PylithScalar scalar = h5.readAttribute<PylithScalar>("/data", "myscalar");
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

    h5.createDataset<int>("/", "data", HDF5::DatasetShape{{{3, 2}}, {{1, 2}}});
    h5.close();

    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);

    h5.open("test.h5", H5F_ACC_RDONLY);
#if defined(PYLITH_HDF5_USE_API_18)
    group.id = H5Gopen2(h5._file, "/", H5P_DEFAULT);
#else
    group.id = H5Gopen(h5._file, "/");
#endif
    REQUIRE(group.id >= 0);
#if defined(PYLITH_HDF5_USE_API_18)
    dataset.id = H5Dopen2(group.id, "data", H5P_DEFAULT);
#else
    dataset.id = H5Dopen(group.id, "data");
#endif
    h5.close();

    CHECK(dataset.id >= 0);

    PYLITH_METHOD_END;
} // testCreateDataset


// ------------------------------------------------------------------------------------------------
// Test writeDatasetChunk() and readDatasetChunk().
void
pylith::meshio::TestHDF5::testDatasetChunk(void) {
    PYLITH_METHOD_BEGIN;

    const int ndimsE = 3;
    const hsize_t dimsE[ndimsE] = { 4, 2, 3 };

    // Create data.
    hsize_t nitems = dimsE[0];
    hsize_t nitemsS = 1;
    for (int i = 1; i < ndimsE; ++i) {
        nitems *= dimsE[i];
        nitemsS *= dimsE[i];
    } // for
    std::unique_ptr<int[]> valuesE((nitems > 0) ? new int[nitems] : nullptr);
    for (size_t i = 0; i < nitems; ++i) {
        valuesE[i] = 2 * i + 1;
    } // for

    HDF5 h5("test.h5", H5F_ACC_TRUNC);
    h5.createDataset<int>("/", "data", HDF5::DatasetShape{{{4, 2, 3}}, {{1, 2, 3}}});

    for (size_t i = 0; i < dimsE[0]; ++i) {
        const std::vector<int> chunk(valuesE.get() + i*nitemsS, valuesE.get() + (i+1)*nitemsS);
        HDF5::ChunkInfo chunkInfo{{i+1, 2, 3}, {1, 2, 3}, (int)i};
        h5.writeDatasetChunk<int>("/", "data", chunk, chunkInfo);
    } // for
    h5.close();

    h5.open("test.h5", H5F_ACC_RDONLY);
    for (size_t i = 0; i < dimsE[0]; ++i) {
        HDF5::Dataset<int> chunk = h5.readDatasetChunk<int>("/", "data", i);
        REQUIRE(ndimsE == chunk.dims.size());
        for (int iDim = 1; iDim < ndimsE; ++iDim) {
            REQUIRE(dimsE[iDim] == chunk.dims[iDim]);
        } // for

        for (size_t ii = 0; ii < nitemsS; ++ii) {
            CHECK(valuesE[i*nitemsS+ii] == chunk.data[ii]);
        }
    } // for

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

    // Create raw data file
    hsize_t nitems = dimsE[0];
    for (int i = 1; i < ndimsE; ++i) {
        nitems *= dimsE[i];
    }
    const hsize_t sizeBytes = nitems * H5Tget_size(H5T_NATIVE_INT);
    std::unique_ptr<int[]> valuesE((nitems > 0) ? new int[nitems] : nullptr);
    for (size_t i = 0; i < nitems; ++i) {
        valuesE[i] = 2 * i + 1;
    } // for
    std::ofstream fout("test.dat");
    fout.write((char*)valuesE.get(), sizeBytes);
    fout.close();

    HDF5 h5("test.h5", H5F_ACC_TRUNC);
    const hid_t intType = sizeof(int) == 64 ? H5T_NATIVE_INT64 : H5T_NATIVE_INT32;
    h5.createDatasetRawExternal("/", "data", "test.dat", HDF5::DatasetDims{{{H5S_UNLIMITED, dimsE[1]}}}, intType);
    h5.extendDatasetRawExternal("/", "data", HDF5::DatasetDims{{{dimsE[0], dimsE[1]}}});
    h5.close();

    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);
    _HDF5::H5Handle dataspace(-1, H5Sclose);

    h5.open("test.h5", H5F_ACC_RDONLY);
#if defined(PYLITH_HDF5_USE_API_18)
    group.id = H5Gopen2(h5._file, "/", H5P_DEFAULT);
#else
    group.id = H5Gopen(h5._file, "/");
#endif
    CHECK(group.id >= 0);
#if defined(PYLITH_HDF5_USE_API_18)
    dataset.id = H5Dopen2(group.id, "data", H5P_DEFAULT);
#else
    dataset.id = H5Dopen(group.id, "data");
#endif
    CHECK(dataset.id >= 0);

    dataspace.id = H5Dget_space(dataset.id);
    CHECK(dataspace.id >= 0);

    const int ndims = H5Sget_simple_extent_ndims(dataspace.id);
    REQUIRE(ndimsE == ndims);
    hsize_t dims[ndimsE];
    H5Sget_simple_extent_dims(dataspace.id, dims, 0);
    for (int i = 0; i < ndims; ++i) {
        CHECK(dimsE[i] == dims[i]);
    } // for

    std::unique_ptr<int[]> values((nitems > 0) ? new int[nitems] : nullptr);
    herr_t err = H5Dread(dataset.id, H5T_NATIVE_INT, dataspace.id, dataspace.id, H5P_DEFAULT, (void*)values.get());
    h5.close();
    CHECK(err >= 0);

    for (size_t i = 0; i < nitems; ++i) {
        CHECK(valuesE[i] == values[i]);
    } // for

    PYLITH_METHOD_END;
} // testDatasetRawExternal


// ------------------------------------------------------------------------------------------------
// Test writeAttribute(string) and readAttribute(string).
void
pylith::meshio::TestHDF5::testAttributeString(void) {
    PYLITH_METHOD_BEGIN;

    HDF5 h5("test.h5", H5F_ACC_TRUNC);

    h5.createDataset<int>("/", "data", HDF5::DatasetShape{{2}, {2}});

    const std::string valueE = "abcd";
    h5.writeAttributeString("/data", "mystring", valueE.c_str());
    h5.close();

    h5.open("test.h5", H5F_ACC_RDONLY);
    std::string value = h5.readAttributeString("/data", "mystring");
    h5.close();
    CHECK(valueE == value);

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
