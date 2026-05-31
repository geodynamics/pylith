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
#include "pylith/utils/journals.hh" // USES journal macros
#include "pylith/utils/Exceptions.hh" // USES Exception

#include <cstring> // USES strlen(), strnlen(), strncpy()
#include <cassert> // USES assert()


// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::meshio::HDF5::HDF5(void) :
    _file(-1) {}


// ------------------------------------------------------------------------------------------------
// Constructor with filename and mode.
pylith::meshio::HDF5::HDF5(const char* filename,
                           hid_t mode) :
    _file(-1) {
    PYLITH_METHOD_BEGIN;

    open(filename, mode);

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::meshio::HDF5::~HDF5(void) noexcept {
    try {
        close();
    } catch (...) {}
}


// ------------------------------------------------------------------------------------------------
// Open HDF5 file.
void
pylith::meshio::HDF5::open(const char* filename,
                           hid_t mode) {
    PYLITH_METHOD_BEGIN;

    assert(filename);

    if (_file >= 0) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output, "HDF5 file already open.");
    }

    if (hid_t(H5F_ACC_TRUNC) == mode) {
        _file = H5Fcreate(filename, mode, H5P_DEFAULT, H5P_DEFAULT);
        if (_file < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                         "Could not create HDF5 file '" << filename << "'.");
        }
    } else {
        _file = H5Fopen(filename, mode, H5P_DEFAULT);
        if (_file < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                         "Could not open existing HDF5 file '" << filename << "'.");
        }
    }

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Close HDF5 file.
void
pylith::meshio::HDF5::close(void) {
    PYLITH_METHOD_BEGIN;

    if (_file >= 0) {
        if (H5Fclose(_file) < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not close HDF5 file.");
        }
    }
    _file = -1;

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Check if HDF5 file is open.
bool
pylith::meshio::HDF5::isOpen(void) const {
    return _file >= 0;
}


// ------------------------------------------------------------------------------------------------
// Check if HDF5 file has a group.
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
        exists = (H5O_TYPE_GROUP == info.type);
        err = H5Oclose(obj);
        assert(err >= 0);
    }

    PYLITH_METHOD_RETURN(exists);
}


// ------------------------------------------------------------------------------------------------
// Check if HDF5 file has a dataset.
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
        exists = (H5O_TYPE_DATASET == info.type);
        err = H5Oclose(obj);
        assert(err >= 0);
    }

    PYLITH_METHOD_RETURN(exists);
}


// ------------------------------------------------------------------------------------------------
// Get dataset dimensions.
pylith::meshio::HDF5::DatasetDims
pylith::meshio::HDF5::getDatasetDims(const char* parent,
                                     const char* name) {
    PYLITH_METHOD_BEGIN;

    assert(isOpen());
    assert(parent);
    assert(name);

    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);
    _HDF5::H5Handle dataspace(-1, H5Sclose);

    DatasetDims result;
    try {
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(_file, parent);
#endif
        if (group.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not open group.");
        }

#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(group.id, name, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(group.id, name);
#endif
        if (dataset.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not open dataset.");
        }

        dataspace.id = H5Dget_space(dataset.id);
        if (dataspace.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not get dataspace.");
        }

        const int ndims = H5Sget_simple_extent_ndims(dataspace.id);
        if (ndims < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not get number of dataset dimensions.");
        }

        result.dims.resize(ndims);
        if (ndims > 0) {
            if (H5Sget_simple_extent_dims(dataspace.id, result.dims.data(), nullptr) < 0) {
                PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not get dataset dimensions.");
            }
        }

    } catch (pylith::Error& err) {
        err.addContext(pylith::ErrorMessage() << "Error occurred while reading dataset '"
                                              << parent << "/" << name << "'.\n");
        throw;
    } catch (const std::exception& err) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Error occurred while reading dataset '"
                     << parent << "/" << name << "':\n" << err.what());
    } catch (...) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Unknown error occurred while reading dataset '"
                     << parent << "/" << name << "'.");
    }

    PYLITH_METHOD_RETURN(result);
}


// ------------------------------------------------------------------------------------------------
// Get names of datasets in a group.
pylith::string_vector
pylith::meshio::HDF5::getGroupDatasets(const char* parent) {
    PYLITH_METHOD_BEGIN;

    assert(isOpen());
    assert(parent);

    _HDF5::H5Handle group(-1, H5Gclose);

    string_vector names;
    try {
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(_file, parent);
#endif
        if (group.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not open group.");
        }

        H5G_info_t ginfo;
        if (H5Gget_info(group.id, &ginfo) < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not get group info.");
        }

        const int gsize = static_cast<int>(ginfo.nlinks);
        names.resize(gsize);
        for (int i = 0; i < gsize; ++i) {
            char buffer[256];
            ssize_t namelen = H5Lget_name_by_idx(group.id, ".", H5_INDEX_NAME,
                                                 H5_ITER_NATIVE, i, buffer,
                                                 sizeof(buffer), H5P_DEFAULT);
            assert(namelen > 0);
            names[i] = buffer;
        }

    } catch (pylith::Error& err) {
        err.addContext(pylith::ErrorMessage() << "Error occurred while getting names of datasets for group '"
                                              << parent << "'.\n");
        throw;
    } catch (const std::exception& err) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Error occurred while getting names of datasets for group '"
                     << parent << "':\n" << err.what());
    } catch (...) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Unknown error occurred while getting names of datasets for group '"
                     << parent << "'.");
    }

    PYLITH_METHOD_RETURN(names);
}


// ------------------------------------------------------------------------------------------------
// Create group.
void
pylith::meshio::HDF5::createGroup(const char* name) {
    PYLITH_METHOD_BEGIN;

    assert(name);

    _HDF5::H5Handle group(-1, H5Gclose);

#if defined(PYLITH_HDF5_USE_API_18)
    group.id = H5Gcreate2(_file, name, 0, H5P_DEFAULT, H5P_DEFAULT);
#else
    group.id = H5Gcreate(_file, name, 0);
#endif
    if (group.id < 0) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Could not create group '" << name << "'.");
    }

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Member: write string attribute to _file.
void
pylith::meshio::HDF5::writeAttributeString(const char* parent,
                                           const char* name,
                                           const char* value) {
    PYLITH_METHOD_BEGIN;

    HDF5::writeAttributeString(_file, parent, name, value);

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Static: write string attribute to an external HDF5 handle.
void
pylith::meshio::HDF5::writeAttributeString(hid_t h5,
                                           const char* parent,
                                           const char* name,
                                           const char* value) {
    PYLITH_METHOD_BEGIN;

    assert(h5 >= 0);
    assert(parent);
    assert(name);
    assert(value);

    _HDF5::H5Handle dataset(-1, H5Dclose);
    _HDF5::H5Handle dataspace(-1, H5Sclose);
    _HDF5::H5Handle datatype(-1, H5Tclose);
    _HDF5::H5Handle attribute(-1, H5Aclose);

    try {
#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(h5, parent, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(h5, parent);
#endif
        if (dataset.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not open parent dataset for");
        }

        dataspace.id = H5Screate(H5S_SCALAR);
        if (dataspace.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not create dataspace for");
        }

        datatype.id = H5Tcopy(H5T_C_S1);
        if (datatype.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not create datatype for");
        }

        if (H5Tset_size(datatype.id, strlen(value) + 1) < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not set size of");
        }

#if defined(PYLITH_HDF5_USE_API_18)
        attribute.id = H5Acreate2(dataset.id, name, datatype.id, dataspace.id,
                                  H5P_DEFAULT, H5P_DEFAULT);
#else
        attribute.id = H5Acreate(dataset.id, name, datatype.id, dataspace.id, H5P_DEFAULT);
#endif
        if (attribute.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not create");
        }

        if (H5Awrite(attribute.id, datatype.id, value) < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not write");
        }

    } catch (pylith::Error& err) {
        err.addContext(pylith::ErrorMessage() << " attribute '" << name << "' of '" << parent << "'.");
        throw;
    } catch (const std::exception& err) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     err.what() << " attribute '" << name << "' of '" << parent << "'.");
    }

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Read string attribute.
std::string
pylith::meshio::HDF5::readAttributeString(const char* parent,
                                          const char* name) {
    PYLITH_METHOD_BEGIN;

    assert(parent);
    assert(name);

    _HDF5::H5Handle dataset(-1, H5Dclose);
    _HDF5::H5Handle attribute(-1, H5Aclose);
    _HDF5::H5Handle datatype(-1, H5Tclose);

    std::string value;
    try {
#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(_file, parent, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(_file, parent);
#endif
        if (dataset.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not open parent dataset for");
        }

        attribute.id = H5Aopen(dataset.id, name, H5P_DEFAULT);
        if (attribute.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not open");
        }

        datatype.id = H5Aget_type(attribute.id);
        if (datatype.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not get datatype of");
        }

        const int len = static_cast<int>(H5Tget_size(datatype.id));
        if (len <= 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Nonpositive size for datatype of");
        }

        std::vector<char> buffer(len);
        if (H5Aread(attribute.id, datatype.id, buffer.data()) < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not read");
        }
        value.assign(buffer.data(), strnlen(buffer.data(), static_cast<size_t>(len)));

    } catch (pylith::Error& err) {
        err.addContext(pylith::ErrorMessage() << " attribute '" << name << "' of '" << parent << "'.");
        throw;
    } catch (const std::exception& err) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     err.what() << " attribute '" << name << "' of '" << parent << "'.");
    }

    PYLITH_METHOD_RETURN(value);
}


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::HDF5::createDatasetRawExternal(const char* parent,
                                               const char* name,
                                               const char* filename,
                                               const DatasetDims& maxDims,
                                               const hid_t datatype) {
    PYLITH_METHOD_BEGIN;

    assert(parent);
    assert(name);
    assert(filename);
    assert(!maxDims.dims.empty());

    const int ndims = static_cast<int>(maxDims.dims.size());

    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle dataspace(-1, H5Sclose);
    _HDF5::H5Handle property(-1, H5Pclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);

    try {
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(_file, parent);
#endif
        if (group.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not open group.");
        }

        std::vector<hsize_t> curDims(maxDims.dims);
        if (curDims[0] == H5S_UNLIMITED) { curDims[0] = 1; }

        dataspace.id = H5Screate_simple(ndims, curDims.data(), maxDims.dims.data());
        if (dataspace.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not create dataspace.");
        }

        property.id = H5Pcreate(H5P_DATASET_CREATE);
        if (property.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not create property for dataset.");
        }

        if (H5Pset_external(property.id, filename, 0, H5F_UNLIMITED) < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not set external file property.");
        }

#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dcreate2(group.id, name, datatype, dataspace.id,
                                H5P_DEFAULT, property.id, H5P_DEFAULT);
#else
        dataset.id = H5Dcreate(group.id, name, datatype, dataspace.id, property.id);
#endif
        if (dataset.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not create dataset.");
        }

    } catch (pylith::Error& err) {
        err.addContext(pylith::ErrorMessage() << "Error occurred while creating dataset '"
                                              << parent << "/" << name << "'.\n");
        throw;
    } catch (const std::exception& err) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Error occurred while creating dataset '"
                     << parent << "/" << name << "':\n" << err.what());
    } catch (...) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Unknown error occurred while creating dataset '" << name << "'.");
    }

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::HDF5::extendDatasetRawExternal(const char* parent,
                                               const char* name,
                                               const DatasetDims& dims) {
    PYLITH_METHOD_BEGIN;

    assert(parent);
    assert(name);
    assert(!dims.dims.empty());

    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);

    try {
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(_file, parent);
#endif
        if (group.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not open group.");
        }

#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(group.id, name, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(group.id, name);
#endif
        if (dataset.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not open dataset.");
        }

#if defined(PYLITH_HDF5_USE_API_18)
        herr_t err = H5Dset_extent(dataset.id, dims.dims.data());
#else
        herr_t err = H5Dextend(dataset.id, dims.dims.data());
#endif
        if (err < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not set dataset extent.");
        }

    } catch (pylith::Error& err) {
        err.addContext(pylith::ErrorMessage() << "Error occurred while updating dataset '"
                                              << parent << "/" << name << "'.\n");
        throw;
    } catch (const std::exception& err) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Error occurred while updating dataset '"
                     << parent << "/" << name << "':\n" << err.what());
    } catch (...) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Unknown error occurred while updating dataset '" << name << "'.");
    }

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Write dataset of strings — packs into fixed-length rows then calls the static.
void
pylith::meshio::HDF5::writeDataset(const char* parent,
                                   const char* name,
                                   const char* const* sarray,
                                   int nstrings) {
    PYLITH_METHOD_BEGIN;

    assert(parent);
    assert(name);
    assert(sarray);
    assert(nstrings > 0);

    int slen = 0;
    for (int i = 0; i < nstrings; ++i) {
        slen = std::max(slen, static_cast<int>(strlen(sarray[i])));
    }
    slen += 1; // null terminator

    try {
        std::vector<char> strfixedlen(nstrings * slen, '\0');
        for (int i = 0; i < nstrings; ++i) {
            strncpy(strfixedlen.data() + i * slen, sarray[i], slen - 1);
        }

        HDF5::writeDataset(_file, parent, name, strfixedlen.data(), nstrings, slen);

    } catch (pylith::Error& err) {
        err.addContext(pylith::ErrorMessage() << "Error occurred while writing dataset '"
                                              << parent << "/" << name << "'.\n");
        throw;
    } catch (const std::exception& err) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Error occurred while writing dataset '"
                     << parent << "/" << name << "':\n" << err.what());
    } catch (...) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Unknown error occurred while writing dataset '" << name << "'.");
    }

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Write dataset comprised of an array of fixed length strings.
void
pylith::meshio::HDF5::writeDataset(const char* parent,
                                   const char* name,
                                   const char* sarray,
                                   int nstrings,
                                   int slen) {
    HDF5::writeDataset(_file, parent, name, sarray, nstrings, slen);
}


// ------------------------------------------------------------------------------------------------
// Static: write fixed-length string dataset to an external HDF5 handle.
void
pylith::meshio::HDF5::writeDataset(hid_t h5,
                                   const char* parent,
                                   const char* name,
                                   const char* sarray,
                                   int nstrings,
                                   int slen) {
    PYLITH_METHOD_BEGIN;

    assert(h5 >= 0);
    assert(parent);
    assert(name);
    assert(sarray);
    assert(nstrings > 0);
    assert(slen > 0);

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
        if (group.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not open group.");
        }

        const hsize_t dims[1] = {static_cast<hsize_t>(nstrings)};
        dataspace.id = H5Screate_simple(1, dims, nullptr);
        if (dataspace.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not create dataspace.");
        }

        datatype.id = H5Tcopy(H5T_C_S1);
        if (datatype.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not create datatype.");
        }
        if (H5Tset_size(datatype.id, static_cast<size_t>(slen)) < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not set size of datatype.");
        }

#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dcreate2(group.id, name, datatype.id, dataspace.id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
        dataset.id = H5Dcreate(group.id, name, datatype.id, dataspace.id, H5P_DEFAULT);
#endif
        if (dataset.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not create dataset.");
        }

        if (H5Dwrite(dataset.id, datatype.id, H5S_ALL, H5S_ALL, H5P_DEFAULT, sarray) < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not write dataset.");
        }

    } catch (pylith::Error& err) {
        err.addContext(pylith::ErrorMessage() << "Error occurred while creating dataset '"
                                              << parent << "/" << name << "'.\n");
        throw;
    } catch (const std::exception& err) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Error occurred while creating dataset '"
                     << parent << "/" << name << "':\n" << err.what());
    } catch (...) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Unknown error occurred while writing dataset '" << name << "'.");
    }

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Read dataset of strings.
pylith::string_vector
pylith::meshio::HDF5::readDataset(const char* parent,
                                  const char* name) {
    PYLITH_METHOD_BEGIN;

    assert(_file >= 0);
    assert(parent);
    assert(name);

    _HDF5::H5Handle group(-1, H5Gclose);
    _HDF5::H5Handle dataset(-1, H5Dclose);
    _HDF5::H5Handle datatype(-1, H5Tclose);
    _HDF5::H5Handle dataspace(-1, H5Sclose);

    string_vector data;
    try {
#if defined(PYLITH_HDF5_USE_API_18)
        group.id = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
        group.id = H5Gopen(_file, parent);
#endif
        if (group.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not open group.");
        }

#if defined(PYLITH_HDF5_USE_API_18)
        dataset.id = H5Dopen2(group.id, name, H5P_DEFAULT);
#else
        dataset.id = H5Dopen(group.id, name);
#endif
        if (dataset.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not open dataset.");
        }

        datatype.id = H5Dget_type(dataset.id);
        if (datatype.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not get datatype.");
        }
        const int slen = static_cast<int>(H5Tget_size(datatype.id));
        if (slen <= 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not get size of datatype.");
        }

        dataspace.id = H5Dget_space(dataset.id);
        if (dataspace.id < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not get dataspace.");
        }
        if (H5Sget_simple_extent_ndims(dataspace.id) != 1) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Expected 1 dimension for string dataset.");
        }
        hsize_t dims[1];
        H5Sget_simple_extent_dims(dataspace.id, dims, nullptr);
        const int nstrings = static_cast<int>(dims[0]);
        if (nstrings <= 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Zero size for dataset.");
        }

        std::vector<char> strfixedlen(nstrings * slen);
        if (H5Dread(dataset.id, datatype.id, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, strfixedlen.data()) < 0) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::output, "Could not read dataset.");
        }

        data.resize(nstrings);
        for (int i = 0; i < nstrings; ++i) {
            const char* p = strfixedlen.data() + i * slen;
            data[i] = std::string(p, strnlen(p, static_cast<size_t>(slen)));
        }

    } catch (pylith::Error& err) {
        err.addContext(pylith::ErrorMessage() << "Error occurred while reading dataset '"
                                              << parent << "/" << name << "'.\n");
        throw;
    } catch (const std::exception& err) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Error occurred while reading dataset '"
                     << parent << "/" << name << "':\n" << err.what());
    } catch (...) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::output,
                     "Unknown error occurred while reading dataset '" << name << "'.");
    }

    PYLITH_METHOD_RETURN(data);
}


// End of file
