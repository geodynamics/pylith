// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/meshio/HDF5.hh" // implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include <cstring> // USES strlen()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

#if H5_VERSION_GE(1,12,0)
#define PYLITH_HDF5_USE_API_112
#endif

#if H5_VERSION_GE(1,8,0)
#define PYLITH_HDF5_USE_API_18
#endif

namespace pylith {
    namespace meshio {
        namespace _HDF5 {
            class H5Handle {
                using closer_t = herr_t (*)(hid_t);
public:

                hid_t id{-1};
                closer_t closer{nullptr};

                H5Handle() = default;
                H5Handle(hid_t h,
                         closer_t c) : id(h), closer(c) {}


                ~H5Handle() {
                    if ((id >= 0) && closer) { closer(id);} }

                H5Handle(const H5Handle&) = delete;
                H5Handle& operator=(const H5Handle&) = delete;

                H5Handle(H5Handle&& o) noexcept : id(o.id), closer(o.closer) {
                    o.id = -1;
                }

                H5Handle& operator=(H5Handle&& o) noexcept {
                    if (this != &o) {
                        if ((id >= 0) && closer) {
                            closer(id);
                        }
                        id = o.id;closer = o.closer;o.id = -1;
                    }
                    return *this;
                }

            };
        }
    }
}

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::meshio::HDF5::HDF5(void) :
    _file(-1) { // constructor
} // constructor


// ------------------------------------------------------------------------------------------------
// Constructor with filename and mode.
pylith::meshio::HDF5::HDF5(const char* filename,
                           hid_t mode) {
    PYLITH_METHOD_BEGIN;

    if (hid_t(H5F_ACC_TRUNC) == mode) {
        _file = H5Fcreate(filename, mode, H5P_DEFAULT, H5P_DEFAULT);
        if (_file < 0) {
            std::ostringstream msg;
            msg << "Could not create HDF5 file '" << filename << "'.";
            throw std::runtime_error(msg.str());
        } // if

    } else {
        _file = H5Fopen(filename, mode, H5P_DEFAULT);
        if (_file < 0) {
            std::ostringstream msg;
            msg << "Could not open existing HDF5 file '" << filename << "'.";
            throw std::runtime_error(msg.str());
        } // if
    } // if/else

    PYLITH_METHOD_END;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::HDF5::~HDF5(void) noexcept {
    try {
        close();
    } catch (...) {}
} // destructor


// ------------------------------------------------------------------------------------------------
// Open HDF5 file.
void
pylith::meshio::HDF5::open(const char* filename,
                           hid_t mode) {
    PYLITH_METHOD_BEGIN;

    assert(filename);

    if (_file >= 0) {
        throw std::runtime_error("HDF5 file already open.");
    } // if

    if (hid_t(H5F_ACC_TRUNC) == mode) {
        _file = H5Fcreate(filename, mode, H5P_DEFAULT, H5P_DEFAULT);
        if (_file < 0) {
            std::ostringstream msg;
            msg << "Could not create HDF5 file '" << filename << "'.";
            throw std::runtime_error(msg.str());
        } // if

    } else {
        _file = H5Fopen(filename, mode, H5P_DEFAULT);
        if (_file < 0) {
            std::ostringstream msg;
            msg << "Could not open existing HDF5 file '" << filename << "'.";
            throw std::runtime_error(msg.str());
        } // if
    } // if/else

    PYLITH_METHOD_END;
} // constructor


// ------------------------------------------------------------------------------------------------
// Close HDF5 file.
void
pylith::meshio::HDF5::close(void) {
    PYLITH_METHOD_BEGIN;

    if (_file >= 0) {
        herr_t err = H5Fclose(_file);
        if (err < 0) {
            throw std::runtime_error("Could not close HDF5 file.");
        }
    } // if
    _file = -1;

    PYLITH_METHOD_END;
} // close


// ------------------------------------------------------------------------------------------------
// Check if HDF5 file is open.
bool
pylith::meshio::HDF5::isOpen(void) const {
    return (_file == -1) ? false : true;
} // isOpen


// ------------------------------------------------------------------------------------------------
// Check if HDF5 file has group.
bool
pylith::meshio::HDF5::hasGroup(const char* name) {
    PYLITH_METHOD_BEGIN;

    assert(isOpen());
    assert(name);

    bool exists = false;
    if (H5Lexists(_file, name, H5P_DEFAULT)) {
        hid_t obj = H5Oopen(_file, name, H5P_DEFAULT);
        assert(obj >= 0);
        H5O_info_t info;
#if defined(PYLITH_HDF5_USE_API_112)
        herr_t err = H5Oget_info(obj, &info, H5O_INFO_ALL);
#else
        herr_t err = H5Oget_info(obj, &info);
#endif
        assert(err >= 0);
        if (H5O_TYPE_GROUP == info.type) {
            exists = true;
        }
        err = H5Oclose(obj);
        assert(err >= 0);
    } // if

    PYLITH_METHOD_RETURN(exists);
} // hasGroup


// ------------------------------------------------------------------------------------------------
// Check if HDF5 file has dataset.
bool
pylith::meshio::HDF5::hasDataset(const char* name) {
    PYLITH_METHOD_BEGIN;

    assert(isOpen());
    assert(name);

    bool exists = false;
    if (H5Lexists(_file, name, H5P_DEFAULT)) {
        hid_t obj = H5Oopen(_file, name, H5P_DEFAULT);
        assert(obj >= 0);
        H5O_info_t info;
#if defined(PYLITH_HDF5_USE_API_112)
        herr_t err = H5Oget_info(obj, &info, H5O_INFO_ALL);
#else
        herr_t err = H5Oget_info(obj, &info);
#endif
        assert(err >= 0);
        if (H5O_TYPE_DATASET == info.type) {
            exists = true;
        }
        err = H5Oclose(obj);
        assert(err >= 0);
    } // if

    PYLITH_METHOD_RETURN(exists);
} // hasDataset


// ------------------------------------------------------------------------------------------------
// Get topology metadata.
void
pylith::meshio::HDF5::getDatasetDims(hsize_t** dims,
                                     int* ndims,
                                     const char* parent,
                                     const char* name) {
    PYLITH_METHOD_BEGIN;

    assert(dims);
    assert(ndims);
    assert(isOpen());

    std::unique_ptr<hsize_t[]> tmpDims;

    // Scope guards to ensure cleanup.
    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);
    _HDF5::H5Handle dataspace(-1, H5Sclose);

    try {
        // Open group
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(_file, parent);
#endif
        if (group.id < 0) {
            throw std::runtime_error("Could not open group.");
        }

        // Open the dataset
#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(group.id, name, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(group, name);
#endif
        if (dataset.id < 0) {
            throw std::runtime_error("Could not open dataset.");
        }

        // Get dataspace
        dataspace.id = H5Dget_space(dataset.id);
        if (dataspace.id < 0) {
            throw std::runtime_error("Could not get dataspace.");
        }

        const int tmpNumDims = H5Sget_simple_extent_ndims(dataspace.id);
        if (tmpNumDims < 0) {
            throw std::runtime_error("Could not get number of dataset dimensions.");
        }
        if (tmpNumDims > 0) {
            tmpDims.reset(new hsize_t[tmpNumDims]);
            if (H5Sget_simple_extent_dims(dataspace.id, tmpDims.get(), nullptr) < 0) {
                throw std::runtime_error("Could not get dataset dimensions.");
            } // if
        } else {
            // Scalar (tmpDims == nullptr)
        } // if/else

        *ndims = tmpNumDims;
        if (*dims) { delete[] *dims; }
        *dims = tmpDims.release();
    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while reading dataset '"
            << parent << "/" << name << "':\n"
            << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        std::ostringstream msg;
        msg << "Unknown error occurred while reading dataset '"
            << parent << "/" << name << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // getDatasetDims


// ------------------------------------------------------------------------------------------------
// Get names of datasets in group.
void
pylith::meshio::HDF5::getGroupDatasets(string_vector* names,
                                       const char* parent) {
    PYLITH_METHOD_BEGIN;

    assert(names);
    assert(isOpen());

    // Scope guards to ensure cleanup.
    _HDF5::H5Handle group(-1, H5Gclose);

    try {
        // Open group
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(_file, parent);
#endif
        if (group.id < 0) {
            throw std::runtime_error("Could not open group.");
        }

        H5G_info_t ginfo;
        herr_t err = H5Gget_info(group.id, &ginfo);
        if (err < 0) {
            throw std::runtime_error("Could not get group info.");
        }
        const int gsize = ginfo.nlinks;

        names->resize(gsize);
        for (int i = 0, index = 0; i < gsize; ++i) {
            char buffer[256];
            ssize_t namelen = H5Lget_name_by_idx(group.id, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, buffer, 256, H5P_DEFAULT);assert(namelen > 0);
            (*names)[index++] = buffer;
        } // for

    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while getting names of datasets for group '"
            << parent << "':\n"
            << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        std::ostringstream msg;
        msg << "Unknown error occurred while getting names of datasets for group '"
            << parent << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // getGroupDatasets


// ------------------------------------------------------------------------------------------------
// Create group.
void
pylith::meshio::HDF5::createGroup(const char* name) {
    PYLITH_METHOD_BEGIN;

    assert(name);

    _HDF5::H5Handle group(-1, H5Gclose);

#if defined(PYLITH_HDF5_USE_API_18)
    group.id = H5Gcreate2(_file, name, 0, H5P_DEFAULT, H5P_DEFAULT);
#else // deprecated HDF5 1.6 API
    group.id = H5Gcreate(_file, name, 0);
#endif
    if (group.id < 0) {
        std::ostringstream msg;
        msg << "Could not create group '" << name << "'.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // createGroup


// ------------------------------------------------------------------------------------------------
// Write scalar attribute.
void
pylith::meshio::HDF5::writeAttribute(const char* parent,
                                     const char* name,
                                     const void* value,
                                     hid_t datatype) {
    PYLITH_METHOD_BEGIN;

    HDF5::writeAttribute(_file, parent, name, value, datatype);

    PYLITH_METHOD_END;
} // writeAttribute


// ------------------------------------------------------------------------------------------------
// Write scalar attribute (external HDF5 handle).
void
pylith::meshio::HDF5::writeAttribute(hid_t h5,
                                     const char* parent,
                                     const char* name,
                                     const void* value,
                                     hid_t datatype) {
    PYLITH_METHOD_BEGIN;

    assert(parent);
    assert(name);
    assert(value);

    // Scope guards to ensure cleanup.
    _HDF5::H5Handle dataspace(-1, H5Sclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);
    _HDF5::H5Handle attribute(-1, H5Aclose);

    try {
        dataspace.id = H5Screate(H5S_SCALAR);
        if (dataspace.id < 0) {
            throw std::runtime_error("Could not create dataspace for");
        }

#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(h5, parent, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(h5, parent);
#endif
        if (dataset.id < 0) {
            throw std::runtime_error("Could not open parent dataset for");
        }

#if defined(PYLITH_HDF5_USE_API_18)
        attribute.id = H5Acreate2(dataset.id, name, datatype, dataspace.id, H5P_DEFAULT, H5P_DEFAULT);
#else
        attribute.id = H5Acreate(dataset, name, datatype, dataspace.id, H5P_DEFAULT);
#endif
        if (attribute.id < 0) {
            throw std::runtime_error("Could not create");
        }

        hid_t err = H5Awrite(attribute.id, datatype, value);
        if (err < 0) {
            throw std::runtime_error("Could not write");
        }

    } catch (std::exception& err) {
        std::ostringstream msg;
        msg << err.what() << " attribute '" << name << "' of '" << parent << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // writeAttribute


// ------------------------------------------------------------------------------------------------
// Write string attribute.
void
pylith::meshio::HDF5::writeAttribute(const char* parent,
                                     const char* name,
                                     const char* value) {
    PYLITH_METHOD_BEGIN;

    HDF5::writeAttribute(_file, parent, name, value);

    PYLITH_METHOD_END;
} // writeAttribute


// ------------------------------------------------------------------------------------------------
// Write string attribute (external handle to HDF5 file).
void
pylith::meshio::HDF5::writeAttribute(hid_t h5,
                                     const char* parent,
                                     const char* name,
                                     const char* value) {
    PYLITH_METHOD_BEGIN;

    assert(h5 >= 0);
    assert(parent);
    assert(name);
    assert(value);

    // Scope guards to ensure cleanup.
    _HDF5::H5Handle attribute(-1, H5Aclose);
    _HDF5::H5Handle datatype(-1, H5Tclose);
    _HDF5::H5Handle dataspace(-1, H5Sclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);

    try {
#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(h5, parent, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(h5, parent);
#endif
        if (dataset.id < 0) {
            throw std::runtime_error("Could not open parent dataset for");
        }

        dataspace.id = H5Screate(H5S_SCALAR);
        if (dataspace.id < 0) {
            throw std::runtime_error("Could not create dataspace for");
        }

        datatype.id = H5Tcopy(H5T_C_S1);
        if (datatype.id < 0) {
            throw std::runtime_error("Could not create datatype for");
        }

        herr_t err = H5Tset_size(datatype.id, strlen(value)+1);
        if (err < 0) {
            throw std::runtime_error("Could not set size of");
        }

#if defined(PYLITH_HDF5_USE_API_18)
        attribute.id = H5Acreate2(dataset.id, name,
                                  datatype.id, dataspace.id, H5P_DEFAULT, H5P_DEFAULT);
#else
        attribute.id = H5Acreate(dataset.id, name,
                                 datatype.id, dataspace.id, H5P_DEFAULT);
#endif
        if (attribute.id < 0) {
            throw std::runtime_error("Could not create");
        }

        err = H5Awrite(attribute.id, datatype.id, value);
        if (err < 0) {
            throw std::runtime_error("Could not write");
        }

    } catch (std::exception& err) {
        std::ostringstream msg;
        msg << err.what() << " attribute '" << name << "' of '" << parent << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // writeAttribute


// ------------------------------------------------------------------------------------------------
// Read scalar attribute.
void
pylith::meshio::HDF5::readAttribute(const char* parent,
                                    const char* name,
                                    void* value) {
    PYLITH_METHOD_BEGIN;

    assert(parent);
    assert(name);
    assert(value);

    // Scope guards to ensure cleanup.
    _HDF5::H5Handle attribute(-1, H5Aclose);
    _HDF5::H5Handle datatype(-1, H5Tclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);

    try {
#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(_file, parent, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(_file, parent);
#endif
        if (dataset.id < 0) {
            throw std::runtime_error("Could not open parent dataset for");
        }

        attribute.id = H5Aopen(dataset.id, name, H5P_DEFAULT);
        if (attribute.id < 0) {
            throw std::runtime_error("Could not open");
        }

        datatype.id = H5Aget_type(attribute.id);
        if (datatype.id < 0) {
            throw std::runtime_error("Could not get datatype of");
        }

        hid_t err = H5Aread(attribute.id, datatype.id, value);
        if (err < 0) {
            throw std::runtime_error("Could not read");
        }

    } catch (std::exception& err) {
        std::ostringstream msg;
        msg << err.what() << " attribute '" << name << "' of '" << parent << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // readAttribute


// ------------------------------------------------------------------------------------------------
// Read string attribute.
std::string
pylith::meshio::HDF5::readAttribute(const char* parent,
                                    const char* name) {
    PYLITH_METHOD_BEGIN;

    assert(parent);
    assert(name);

    // Scope guards to ensure cleanup.
    _HDF5::H5Handle attribute(-1, H5Aclose);
    _HDF5::H5Handle datatype(-1, H5Tclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);

    std::string value;
    try {
#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(_file, parent, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(_file, parent);
#endif
        if (dataset.id < 0) {
            throw std::runtime_error("Could not open parent dataset for");
        }

        attribute.id = H5Aopen(dataset.id, name, H5P_DEFAULT);
        if (attribute.id < 0) {
            throw std::runtime_error("Could not open");
        }

        datatype.id = H5Aget_type(attribute.id);
        if (datatype.id < 0) {
            throw std::runtime_error("Could not get datatype of");
        }

        // :TODO: Check that datatype is a string

        const int len = H5Tget_size(datatype.id);
        if (len <= 0) {
            throw std::runtime_error("Nonpositive size for datatype of");
        }

        std::unique_ptr<char[]> buffer((len > 0) ? new char[len] : nullptr);
        hid_t err = H5Aread(attribute.id, datatype.id, (void*)buffer.get());
        if (err < 0) {
            throw std::runtime_error("Could not read");
        }
        value.assign(buffer.get(), strnlen(buffer.get(), len));
        buffer.reset(nullptr);

    } catch (std::exception& err) {
        std::ostringstream msg;
        msg << err.what() << " attribute '" << name << "' of '" << parent << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_RETURN(std::string(value));
} // readAttribute


// ------------------------------------------------------------------------------------------------
// Create dataset.
void
pylith::meshio::HDF5::createDataset(const char* parent,
                                    const char* name,
                                    const hsize_t* maxDims,
                                    const hsize_t* dimsChunk,
                                    const int ndims,
                                    hid_t datatype) {
    PYLITH_METHOD_BEGIN;

    assert(parent);
    assert(name);
    assert(maxDims);
    assert(dimsChunk);

    // Scope guards to ensure cleanup.
    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle dataspace(-1, H5Sclose);
    _HDF5::H5Handle property(-1, H5Pclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);

    try {
        // Open group
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(_file, parent);
#endif
        if (group.id < 0) {
            throw std::runtime_error("Could not open group.");
        }

        // Create the dataspace
        std::unique_ptr<hsize_t[]> curDims((ndims > 0) ? new hsize_t[ndims] : nullptr);
        if (ndims > 0) {
            curDims.get()[0] = (maxDims[0] == H5S_UNLIMITED) ? 1 : maxDims[0];
        }
        for (int i = 1; i < ndims; ++i) {
            curDims.get()[i] = maxDims[i];
        }
        dataspace.id = H5Screate_simple(ndims, curDims.get(), maxDims);
        curDims.reset(nullptr);
        if (dataspace.id < 0) {
            throw std::runtime_error("Could not create dataspace.");
        }

        // Create chunked dataset
        property.id = H5Pcreate(H5P_DATASET_CREATE);
        if (property.id < 0) {
            throw std::runtime_error("Could not create property for dataset.");
        }

        herr_t err = H5Pset_chunk(property.id, ndims, dimsChunk);
        if (err < 0) {
            throw std::runtime_error("Could not set chunk.");
        }

        // Set gzip compression level for chunk.
        H5Pset_deflate(property.id, 6);

#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dcreate2(group.id, name, datatype, dataspace.id, H5P_DEFAULT, property.id, H5P_DEFAULT);
#else
        dataset.id = H5Dcreate(group.id, name, datatype, dataspace.id, property.id);
#endif
        if (dataset.id < 0) {
            throw std::runtime_error("Could not create dataset.");
        }

    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while creating dataset '"
            << parent << "/" << name << "':\n"
            << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        std::ostringstream msg;
        msg << "Unknown error occurred while creating dataset '" << name << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // createDataset


// ------------------------------------------------------------------------------------------------
// Read dataset.
void
pylith::meshio::HDF5::readDataset(char** const data,
                                  hsize_t** const dims,
                                  int* const ndims,
                                  hid_t datatype,
                                  const char* parent,
                                  const char* name) {
    PYLITH_METHOD_BEGIN;

    assert(parent);
    assert(name);
    assert(data);
    assert(dims);
    assert(ndims);
    assert(_file >= 0);

    // Temporary buffers. Pass to output args on success.
    std::unique_ptr<hsize_t[]> tmpDims;
    std::unique_ptr<char[]> tmpData;

    // Scope guards to ensure cleanup.
    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle dataspace(-1, H5Sclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);

    try {
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(_file, parent);
#endif
        if (group.id < 0) {
            throw std::runtime_error("Could not open group.");
        }

        // Open the dataset
#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(group.id, name, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(group.id, name);
#endif
        if (dataset.id < 0) {
            throw std::runtime_error("Could not open dataset.");
        }

        dataspace.id = H5Dget_space(dataset.id);
        if (dataspace.id < 0) {
            throw std::runtime_error("Could not get dataspace.");
        }

        const int tmpNumDims = H5Sget_simple_extent_ndims(dataspace.id);
        if (tmpNumDims < 0) {
            throw std::runtime_error("Could not get number of dimensions.");
        }

        if (tmpNumDims > 0) {
            tmpDims.reset(new hsize_t[tmpNumDims]);
            if (H5Sget_simple_extent_dims(dataspace.id, tmpDims.get(), nullptr) < 0) {
                throw std::runtime_error("could not get dimensions.");
            }
        } else {
            // Scalar (tmpDims == nullptr)
        }

        size_t datatypeSize = H5Tget_size(datatype);
        if (datatypeSize == 0) {
            throw std::runtime_error("Datatype size is zero.");
        }
        size_t totalBytes = datatypeSize;
        for (int i = 0; i < tmpNumDims; ++i) {
            const size_t dim = static_cast<size_t>(tmpDims[i]);
            if ((dim != 0) && (totalBytes > std::numeric_limits<size_t>::max() / dim)) {
                throw std::overflow_error("Dataset size exceeds size_t range.");
            }
            totalBytes *= dim;
        }

        if (totalBytes > 0) {
            tmpData.reset(new char[totalBytes]);
            if (H5Dread(dataset.id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)tmpData.get()) < 0) {
                throw std::runtime_error("Could not read data.");
            }
        } else {
            tmpData.reset(nullptr);
        } // if/else

        *ndims = tmpNumDims;
        if (*dims) { delete[] *dims; }
        if (*data) { delete[] *data; }
        *dims = tmpDims.release();
        *data = tmpData.release();

    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while reading dataset '"
            << parent << "/" << name << "':\n"
            << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        std::ostringstream msg;
        msg << "Unknown error occurred while reading dataset '"
            << parent << "/" << name << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // readDataset


// ------------------------------------------------------------------------------------------------
// Append slice to dataset.
void
pylith::meshio::HDF5::writeDatasetChunk(const char* parent,
                                        const char* name,
                                        const void* data,
                                        const hsize_t* dims,
                                        const hsize_t* dimsChunk,
                                        const int ndims,
                                        const int chunk,
                                        hid_t datatype) { // writeDatasetChunk
    PYLITH_METHOD_BEGIN;

    assert(parent);
    assert(name);
    assert(data);
    assert(dims);
    assert(_file >= 0);

    // Scope guards to ensure cleanup.
    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle dataspace(-1, H5Sclose);
    _HDF5::H5Handle chunkspace(-1, H5Sclose);
    _HDF5::H5Handle property(-1, H5Pclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);

    try {
        // Select hyperslab in file
        std::unique_ptr<hsize_t[]> count((ndims > 0) ? new hsize_t[ndims] : nullptr);
        std::unique_ptr<hsize_t[]> stride((ndims > 0) ? new hsize_t[ndims] : nullptr);
        std::unique_ptr<hsize_t[]> offset((ndims > 0) ? new hsize_t[ndims] : nullptr);
        for (int i = 0; i < ndims; ++i) {
            count[i] = 1;
            stride[i] = 1;
            offset[i] = 0;
        } // for
        offset[0] = chunk;

        // Open group
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(_file, parent);
#endif
        if (group.id < 0) {
            throw std::runtime_error("Could not open group.");
        }

        // Open the dataset
#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(group.id, name, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(group.id, name);
#endif
        if (dataset.id < 0) {
            throw std::runtime_error("Could not open dataset.");
        }

#if defined(PYLITH_HDF5_USE_API_18)
        herr_t err = H5Dset_extent(dataset.id, dims);
#else
        herr_t err = H5Dextend(dataset, dims);
#endif
        if (err < 0) {
            throw std::runtime_error("Could not set dataset extent.");
        }

        dataspace.id = H5Dget_space(dataset.id);
        if (dataspace.id < 0) {
            throw std::runtime_error("Could not get dataspace.");
        }

        chunkspace.id = H5Screate_simple(ndims, dimsChunk, 0);
        if (chunkspace.id < 0) {
            throw std::runtime_error("Could not create chunk dataspace.");
        }

        err = H5Sselect_hyperslab(dataspace.id, H5S_SELECT_SET, offset.get(), stride.get(), count.get(), dimsChunk);
        count.reset(nullptr);
        stride.reset(nullptr);
        offset.reset(nullptr);
        if (err < 0) {
            throw std::runtime_error("Could not select hyperslab.");
        }

        err = H5Dwrite(dataset.id, datatype, chunkspace.id, dataspace.id, H5P_DEFAULT, data);
        if (err < 0) {
            throw std::runtime_error("Could not write data.");
        }

    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while writing dataset '"
            << parent << "/" << name << "':\n"
            << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        std::ostringstream msg;
        msg << "Unknown error occurred while writing dataset '"
            << parent << "/" << name << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // writeDatasetChunk


// ------------------------------------------------------------------------------------------------
// Read dataset slice.
void
pylith::meshio::HDF5::readDatasetChunk(const char* parent,
                                       const char* name,
                                       char** const data,
                                       hsize_t** const dimsChunk,
                                       int* const ndims,
                                       const int chunk,
                                       hid_t datatype) {
    PYLITH_METHOD_BEGIN;

    assert(parent);
    assert(name);
    assert(data);
    assert(dimsChunk);
    assert(_file >= 0);

    // Temporary buffers. Pass to output args on success.
    std::unique_ptr<hsize_t[]> tmpDimsChunk;
    std::unique_ptr<char[]> tmpData;

    // Scope guards to ensure cleanup.
    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle dataspace(-1, H5Sclose);
    _HDF5::H5Handle chunkspace(-1, H5Sclose);
    _HDF5::H5Handle property(-1, H5Pclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);

    try {
        // Open group
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(_file, parent);
#endif
        if (group.id < 0) {
            throw std::runtime_error("Could not open group.");
        }

        // Open the dataset
#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(group.id, name, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(group.id, name);
#endif
        if (dataset.id < 0) {
            throw std::runtime_error("Could not open dataset.");
        }

        dataspace.id = H5Dget_space(dataset.id);
        if (dataspace.id < 0) {
            throw std::runtime_error("Could not get dataspace.");
        }

        const int tmpNumDims = H5Sget_simple_extent_ndims(dataspace.id);
        if (tmpNumDims < 0) {
            throw std::runtime_error("Could not get number of dimensions.");
        }

        std::unique_ptr<hsize_t[]> dims((tmpNumDims > 0) ? new hsize_t[tmpNumDims] : nullptr);
        if (tmpNumDims > 0) {
            tmpDimsChunk.reset(new hsize_t[tmpNumDims]);
            if (H5Sget_simple_extent_dims(dataspace.id, dims.get(), nullptr) < 0) {
                throw std::runtime_error("Could not get dimensions.");
            }
        }

        // Select hyperslab in file
        std::unique_ptr<hsize_t[]> count((tmpNumDims > 0) ? new hsize_t[tmpNumDims] : nullptr);
        std::unique_ptr<hsize_t[]> stride((tmpNumDims > 0) ? new hsize_t[tmpNumDims] : nullptr);
        std::unique_ptr<hsize_t[]> offset((tmpNumDims > 0) ? new hsize_t[tmpNumDims] : nullptr);

        for (int i = 0; i < tmpNumDims; ++i) {
            tmpDimsChunk[i] = dims[i];
            count[i] = 1;
            stride[i] = 1;
            offset[i] = 0;
        } // for
        (tmpDimsChunk)[0] = 1;
        offset[0] = chunk;

        chunkspace.id = H5Screate_simple(tmpNumDims, tmpDimsChunk.get(), 0);
        if (chunkspace.id < 0) {
            throw std::runtime_error("Could not create chunk dataspace.");
        }

        herr_t err = H5Sselect_hyperslab(dataspace.id, H5S_SELECT_SET, offset.get(), stride.get(), count.get(), tmpDimsChunk.get());
        if (err < 0) {
            throw std::runtime_error("Could not select hyperslab.");
        }
        count.reset(nullptr);
        stride.reset(nullptr);
        offset.reset(nullptr);

        size_t datatypeSize = H5Tget_size(datatype);
        if (datatypeSize == 0) {
            throw std::runtime_error("Datatype size is zero.");
        }
        size_t totalBytes = datatypeSize;
        for (int i = 0; i < tmpNumDims; ++i) {
            const size_t dim = static_cast<size_t>(tmpDimsChunk[i]);
            if ((dim != 0) && (totalBytes > std::numeric_limits<size_t>::max() / dim)) {
                throw std::overflow_error("Chunk size exceeds size_t range.");
            }
            totalBytes *= dim;
        }

        if (totalBytes > 0) {
            tmpData.reset(new char[totalBytes]);
            err = H5Dread(dataset.id, datatype, chunkspace.id, dataspace.id, H5P_DEFAULT, (void*)tmpData.get());
            if (err < 0) {
                throw std::runtime_error("Could not read data.");
            }
        } else {
            tmpData.reset(nullptr);
        } // if/else

        *ndims = tmpNumDims;
        if (*dimsChunk) { delete[] *dimsChunk; }
        if (*data) { delete[] *data; }
        *dimsChunk = tmpDimsChunk.release();
        *data = tmpData.release();
    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while reading dataset '"
            << parent << "/" << name << "':\n"
            << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        std::ostringstream msg;
        msg << "Unknown error occurred while reading dataset '"
            << parent << "/" << name << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // readDatasetChunk


// ------------------------------------------------------------------------------------------------
// Create dataset associated with data stored in a raw external binary
// file.
void
pylith::meshio::HDF5::createDatasetRawExternal(const char* parent,
                                               const char* name,
                                               const char* filename,
                                               const hsize_t* maxDims,
                                               const int ndims,
                                               hid_t datatype) {
    PYLITH_METHOD_BEGIN;

    assert(parent);
    assert(name);
    assert(filename);
    assert(maxDims);

    // Scope guards to ensure cleanup.
    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle dataspace(-1, H5Sclose);
    _HDF5::H5Handle property(-1, H5Pclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);

    try {
        // Open group
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(_file, parent);
#endif
        if (group.id < 0) {
            throw std::runtime_error("Could not open group.");
        }

        // Create the dataspace
        std::unique_ptr<hsize_t[]> curDims((ndims > 0) ? new hsize_t[ndims] : nullptr);
        if (ndims > 0) {
            curDims[0] = (maxDims[0] == H5S_UNLIMITED) ? 1 : maxDims[0];
        }
        for (int i = 1; i < ndims; ++i) {
            curDims[i] = maxDims[i];
        }
        dataspace.id = H5Screate_simple(ndims, curDims.get(), maxDims);
        if (dataspace.id < 0) {
            throw std::runtime_error("Could not create dataspace.");
        }
        curDims.reset(nullptr);

        // Create property for external dataset
        property.id = H5Pcreate(H5P_DATASET_CREATE);
        if (property.id < 0) {
            throw std::runtime_error("Could not create property for dataset.");
        }

        // Set external file
        const off_t offset = 0;
        herr_t err = H5Pset_external(property.id, filename, offset, H5F_UNLIMITED);
        if (err < 0) {
            throw std::runtime_error("Could not set external file property.");
        }

#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dcreate2(group.id, name, datatype, dataspace.id, H5P_DEFAULT, property.id, H5P_DEFAULT);
#else
        dataset.id = H5Dcreate(group, name, datatype, dataspace, property);
#endif
        if (dataset.id < 0) {
            throw std::runtime_error("Could not create dataset.");
        }

    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while creating dataset '"
            << parent << "/" << name << "':\n"
            << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        std::ostringstream msg;
        msg << "Unknown error occurred while creating dataset '" << name << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // createDatasetRawExternal


// ------------------------------------------------------------------------------------------------
// Create dataset associated with data stored in a raw external binary
// file.
void
pylith::meshio::HDF5::extendDatasetRawExternal(const char* parent,
                                               const char* name,
                                               const hsize_t* dims,
                                               const int ndims) {
    PYLITH_METHOD_BEGIN;

    assert(parent);
    assert(name);
    assert(dims);

    // Scope guards to ensure cleanup.
    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);
    _HDF5::H5Handle property(-1, H5Pclose);
    try {
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(_file, parent);
#endif
        if (group.id < 0) {
            throw std::runtime_error("Could not open group.");
        }

        // Open dataset.
#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(group.id, name, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(group.id, name);
#endif
        if (dataset.id < 0) {
            throw std::runtime_error("Could not open dataset.");
        }

#if defined(PYLITH_HDF5_USE_API_18)
        herr_t err = H5Dset_extent(dataset.id, dims);
#else
        herr_t err = H5Dextend(dataset.id, dims);
#endif
        if (err < 0) {
            throw std::runtime_error("Could not set dataset extent.");
        }

    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while updating dataset '"
            << parent << "/" << name << "':\n"
            << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        std::ostringstream msg;
        msg << "Unknown error occurred while updating dataset '" << name << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // extendDatasetRawExternal


// ------------------------------------------------------------------------------------------------
// Write dataset comprised of an array of strings (external HDF5 handle).
void
pylith::meshio::HDF5::writeDataset(const char* parent,
                                   const char* name,
                                   const char* const* sarray,
                                   const int nstrings) { // writeDataset
    PYLITH_METHOD_BEGIN;

    int slen = 0;
    for (int i = 0; i < nstrings; ++i) {
        slen = std::max(slen, int(strlen(sarray[i])));
    } // for
    slen += 1; // add space for null terminator.

    try {
        std::unique_ptr<char[]> strfixedlen((nstrings*slen > 0) ? new char[nstrings*slen] : nullptr);
        for (int i = 0; i < nstrings; ++i) {
            const int index = i*slen;
            strncpy(&strfixedlen[index], sarray[i], slen-1);
            strfixedlen[index+slen-1] = '\0';
            // Fill rest of string with null.
            for (int j = strlen(sarray[i]); j < slen; ++j) {
                strfixedlen[index+j] = '\0';
            } // for
        } // for

        HDF5::writeDataset(_file, parent, name, strfixedlen.get(), nstrings, slen);
    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while writing dataset '" << parent << "/" << name << "':\n" << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        std::ostringstream msg;
        msg << "Unknown error occurred while writing dataset '" << name << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // writeDataset


// ------------------------------------------------------------------------------------------------
// Write dataset comprised of an array of strings.
void
pylith::meshio::HDF5::writeDataset(const char* parent,
                                   const char* name,
                                   const char* sarray,
                                   const int nstrings,
                                   const int slen) { // writeDataset
    PYLITH_METHOD_BEGIN;

    HDF5::writeDataset(_file, parent, name, sarray, nstrings, slen);

    PYLITH_METHOD_END;
} // writeDataset


// ------------------------------------------------------------------------------------------------
// Write dataset comprised of an array of strings (external HDF5 handle).
void
pylith::meshio::HDF5::writeDataset(hid_t h5,
                                   const char* parent,
                                   const char* name,
                                   const char* sarray,
                                   const int nstrings,
                                   const int slen) { // writeDataset
    PYLITH_METHOD_BEGIN;

    assert(h5 >= 0);
    assert(parent);
    assert(name);
    assert(sarray);
    assert(nstrings > 0);
    assert(slen > 0);

    // Scope guards to ensure cleanup.
    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle datatype(-1, H5Tclose);
    _HDF5::H5Handle dataspace(-1, H5Sclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);
    try {
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(h5, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(h5, parent);
#endif
        if (group.id < 0) { throw std::runtime_error("Could not open group.");}

        // Create the dataspace
        const int ndims = 1;
        hsize_t dims[ndims];dims[0] = nstrings;
        dataspace.id = H5Screate_simple(ndims, dims, NULL);
        if (dataspace.id < 0) { throw std::runtime_error("Could not create dataspace.");}

        datatype.id = H5Tcopy(H5T_C_S1);
        if (datatype.id < 0) { throw std::runtime_error("Could not create datatype.");}
        herr_t err = H5Tset_size(datatype.id, slen);
        if (err < 0) { throw std::runtime_error("Could not set size of datatype.");}

#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dcreate2(group.id, name, datatype.id, dataspace.id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
        dataset.id = H5Dcreate(group, name, datatype, dataspace, H5P_DEFAULT);
#endif
        if (dataset.id < 0) {
            throw std::runtime_error("Could not create dataset.");
        }

        err = H5Dwrite(dataset.id, datatype.id, H5S_ALL, H5S_ALL, H5P_DEFAULT, sarray);
        if (err < 0) { throw std::runtime_error("Could not write dataset.");}

    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while creating dataset '" << parent << "/" << name << "':\n" << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        std::ostringstream msg;
        msg << "Unknown error occurred while writing dataset '" << name << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // writeDataset


// ------------------------------------------------------------------------------------------------
// Read dataset comprised of an array of strings.
pylith::string_vector
pylith::meshio::HDF5::readDataset(const char* parent,
                                  const char* name) {
    PYLITH_METHOD_BEGIN;

    assert(_file >= 0);
    assert(parent);
    assert(name);

    // Scope guards to ensure cleanup.
    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle datatype(-1, H5Tclose);
    _HDF5::H5Handle dataspace(-1, H5Sclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);

    pylith::string_vector data;
    try {
        // Open group
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(_file, parent);
#endif
        if (group.id < 0) {
            throw std::runtime_error("Could not open group.");
        }

        // Open the dataset
#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(group.id, name, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(group, name);
#endif
        if (dataset.id < 0) {
            throw std::runtime_error("Could not open dataset.");
        }

        // Get the datatype
        datatype.id = H5Dget_type(dataset.id);
        if (datatype.id < 0) {
            throw std::runtime_error("Could not get datatype.");
        }
        const int slen = H5Tget_size(datatype.id);
        if (slen < 0) {
            throw std::runtime_error("Could not get size of datatype.");
        }

        // Get the dataspace
        dataspace.id = H5Dget_space(dataset.id);
        if (dataspace.id < 0) {
            throw std::runtime_error("Could not get dataspace.");
        }
        const int ndims = H5Sget_simple_extent_ndims(dataspace.id);
        if (ndims != 1) {
            throw std::runtime_error("Expected 1 dimension for string dataset.");
        }
        hsize_t dims[1];
        H5Sget_simple_extent_dims(dataspace.id, dims, NULL);
        const int nstrings = dims[0];
        if (nstrings <= 0) {
            throw std::runtime_error("Zero size for dataset.");
        }

        // Read the dataset
        std::unique_ptr<char[]> strfixedlen((nstrings*slen > 0) ? new char[nstrings*slen] : nullptr);
        herr_t err = H5Dread(dataset.id, datatype.id, H5S_ALL, H5S_ALL, H5P_DEFAULT, strfixedlen.get());
        if (err < 0) {
            throw std::runtime_error("Could not read dataest.");
        }

        data.resize(nstrings);
        const char* base = strfixedlen.get();
        for (int i = 0; i < nstrings; ++i) {
            const char* p = base + i * slen;
            data[i] = std::string(p, strnlen(p, slen));
        } // for

    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while reading dataset '"
            << parent << "/" << name << "':\n"
            << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        std::ostringstream msg;
        msg << "Unknown error occurred while reading dataset '" << name << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_RETURN(data);
} // readDataset


// End of file
