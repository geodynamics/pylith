// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/meshio/meshiofwd.hh"
#include "pylith/utils/arrayfwd.hh" // USES string_vector

#include "hdf5.h"

#include <string>
#include <vector>

class pylith::meshio::HDF5 {
    friend class TestHDF5; // Unit testing

    // PUBLIC STRUCTS ---------------------------------------------------------------------
public:

    // Structs for packaging method arguments and return values.

    /// Dataset.
    template <typename T>
    struct Dataset {
        std::vector<T>      data;
        std::vector<hsize_t> dims;
    };

    /// Dimensions of a dataset.
    struct DatasetDims {
        std::vector<hsize_t> dims;
    };

    /// Full shape specification for a chunked dataset.
    struct DatasetShape {
        std::vector<hsize_t> maxDims; ///< Maximum extents; H5S_UNLIMITED allowed in [0].
        std::vector<hsize_t> chunkDims; ///< Chunk extents (same rank as maxDims).
    };

    /// Identifies a single chunk within a chunked dataset.
    struct ChunkInfo {
        std::vector<hsize_t> fullDims; ///< Full dataset extent to set via H5Dset_extent.
        std::vector<hsize_t> chunkDims; ///< Shape of the memory buffer (fullDims but [0]=1).
        int chunk{0}; ///< Index along the first (time) axis.
    };


    // PUBLIC METHODS ---------------------------------------------------------------------
public:

    /// Default constructor.
    HDF5(void);

    /** Constructor with filename and mode.
     *
     * @param filename Name of HDF5 file
     * @param mode Mode for HDF5 file
     */
    HDF5(const char* filename,
         hid_t mode);

    /// Destructor
    ~HDF5() noexcept;

    HDF5(const HDF5&) = delete;
    HDF5& operator=(const HDF5&) = delete;

    /** Open HDF5 file.
     *
     * @param filename Name of HDF5 file
     * @param mode Mode for HDF5 file
     */
    void open(const char* filename,
              hid_t mode);

    /// Close HDF5 file.
    void close();

    /** Check if HDF5 file is open.
     *
     * @returns True if HDF5 file is open, false otherwise.
     */
    bool isOpen() const;

    /** Check if HDF5 file has group.
     *
     * @param name Full name of group.
     * @returns True if group exists, false otherwise.
     */
    bool hasGroup(const char* name);

    /** Check if HDF5 file has dataset.
     *
     * @param name Full name of dataset.
     * @returns True if dataset exists, false otherwise.
     */
    bool hasDataset(const char* name);

    /** Create group.
     *
     * Create group and leave group open for further operations.
     *
     * @param name Name of group (with absolute path).
     */
    void createGroup(const char* name);


    /** Get dimensions of dataset.
     *
     * @param parent Full path of parent dataset for attribute.
     * @param name Name of attribute.
     * @returns Dataset dimensions.
     */
    DatasetDims getDatasetDims(const char* parent,
                               const char* name);

    /** Get names of datasets in group.
     *
     * @param group Name of parent.
     * @returns Names of datasets.
     */
    string_vector getGroupDatasets(const char* parent);

    /** Set scalar attribute.
     *
     * @param parent Full path of parent dataset for attribute.
     * @param name Name of attribute.
     * @param value Attribute value.
     */
    template <typename T>
    void writeAttribute(const char* parent,
                        const char* name,
                        T value);

    /** Read scalar attribute.
     *
     * @param parent Full path of parent dataset for attribute.
     * @param name Name of attribute.
     * @returns Attribute value.
     */
    template <typename T>
    T readAttribute(const char* parent,
                    const char* name);

    /** Set string attribute.
     *
     * @param parent Full path of parent dataset for attribute.
     * @param name Name of attribute.
     * @param value String value
     */
    void writeAttributeString(const char* parent,
                              const char* name,
                              const char* value);

    /** Read string attribute.
     *
     * @param parent Full path of parent dataset for attribute.
     * @param name Name of attribute.
     * @returns value String value
     */
    std::string readAttributeString(const char* parent,
                                    const char* name);

    /** Set scalar attribute (used with external handle to HDF5 file,
     * such as PetscHDF5Viewer).
     *
     * @param h5 HDF5 file.
     * @param parent Full path of parent dataset for attribute.
     * @param name Name of attribute.
     * @param value Attribute value.
     */
    template <typename T>
    static void writeAttribute(hid_t h5,
                               const char* parent,
                               const char* name,
                               T value);

    /** Set string attribute (used with external handle to HDF5 file,
     * such as PetscHDF5Viewer).
     *
     * @param h5 HDF5 file.
     * @param parent Full path of parent dataset for attribute.
     * @param name Name of attribute.
     * @param value String value
     */
    static void writeAttributeString(hid_t h5,
                                     const char* parent,
                                     const char* name,
                                     const char* value);

    /** Create dataset.
     *
     * @param parent Full path of parent group for dataset.
     * @param name Name of dataset.
     * @param shape Dataset dimensions.
     */
    template <typename T>
    void createDataset(const char* parent,
                       const char* name,
                       const DatasetShape& shape);

    /** Read dataset.
     *
     * @param parent Full path of parent group for dataset.
     * @param name Name of dataset.
     * @returns Dataset.
     */
    template <typename T>
    Dataset<T> readDataset(const char* parent,
                           const char* name);

    /** Append chunk to dataset.
     *
     * @param parent Full path of parent group for dataset.
     * @param name Name of dataset.
     * @param data Data.
     * @param chunk Chunk information.
     */
    template <typename T>
    void writeDatasetChunk(const char* parent,
                           const char* name,
                           const std::vector<T>& data,
                           const ChunkInfo& chunk);

    /** Read dataset chunk.
     *
     * Currently this method assumes the chunk size (slice along dim=0).
     *
     * @param parent Full path of parent group for dataset.
     * @param name Name of dataset.
     * @param chunk Index of data slice.
     * @returns Dataset chunk.
     */
    template <typename T>
    Dataset<T> readDatasetChunk(const char* parent,
                                const char* name,
                                int chunk);

    /** Create dataset associated with data stored in a raw external
     * binary file.
     *
     * @param parent Full path of parent group for dataset.
     * @param name Name of dataset.
     * @param filename Name of external raw data file.
     * @param maxDims Maximum dimensions.
     * @param datatype Datatype of external data.
     */
    void createDatasetRawExternal(const char* parent,
                                  const char* name,
                                  const char* filename,
                                  const DatasetDims& maxDims,
                                  const hid_t datatype);

    /** Update the properties of a dataset associated with data stored
     * in a raw external binary file.
     *
     * @param parent Full path of parent group for dataset.
     * @param name Name of dataset.
     * @param dims Dimensions of data.
     */
    void extendDatasetRawExternal(const char* parent,
                                  const char* name,
                                  const DatasetDims& dims);

    /** Write dataset comprised of an array of strings.
     *
     * Data is written as fixed length strings matching the maximum
     * string length.
     *
     * @param parent Full path of parent group for dataset.
     * @param name Name of dataset.
     * @param sarray Array of null terminated C strings.
     * @param nstrings Size of array.
     */
    void writeDataset(const char* parent,
                      const char* name,
                      const char* const* sarray,
                      int nstrings);

    /** Write dataset comprised of an array of fixed length strings.
     *
     * @param parent Full path of parent group for dataset.
     * @param name Name of dataset.
     * @param sarray Array of null terminated C strings.
     * @param nstrings Size of array.
     * @param slen Fixed string length.
     */
    void writeDataset(const char* parent,
                      const char* name,
                      const char* sarray,
                      int nstrings,
                      int slen);

    /** Write dataset comprised of an array of fixed length strings.
     *
     * @param parent Full path of parent group for dataset.
     * @param name Name of dataset.
     * @param sarray Array of null terminated C strings.
     * @param nstrings Size of array.
     * @param slen Fixed string length.
     */
    static void writeDataset(hid_t h5,
                             const char* parent,
                             const char* name,
                             const char* sarray,
                             int nstrings,
                             int slen);

    /** Read dataset comprised of an array of strings.
     *
     * @param parent Full path of parent group for dataset.
     * @param name Name of dataset.
     * @returns Array of string.
     */
    pylith::string_vector readDataset(const char* parent,
                                      const char* name);

    // PRIVATE METHODS --------------------------------------------------------------------
private:

    // Map a C++ scalar type to the corresponding native HDF5 type.
    template <typename T> static hid_t _h5NativeType();

    // PRIVATE MEMBERS --------------------------------------------------------------------
private:

    hid_t _file;
}; // HDF5


#include "HDF5.icc"

// End of file
