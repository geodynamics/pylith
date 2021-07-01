// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "HDF5.hh" // implementation of class methods

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

// ----------------------------------------------------------------------
// Default constructor.
pylith::meshio::HDF5::HDF5(void) :
  _file(-1)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Constructor with filename and mode.
pylith::meshio::HDF5::HDF5(const char* filename,
			   hid_t mode)
{ // constructor
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

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::HDF5::~HDF5(void)
{ // destructor
  close();
} // destructor

// ----------------------------------------------------------------------
// Open HDF5 file.
void
pylith::meshio::HDF5::open(const char* filename,
			   hid_t mode)
{ // open
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

// ----------------------------------------------------------------------
// Close HDF5 file.
void
pylith::meshio::HDF5::close(void)
{ // close
  PYLITH_METHOD_BEGIN;

  if (_file >= 0) {
    herr_t err = H5Fclose(_file);
    if (err < 0)
      throw std::runtime_error("Could not close HDF5 file.");
  } // if
  _file = -1;

  PYLITH_METHOD_END;
} // close

// ----------------------------------------------------------------------
// Check if HDF5 file is open.
bool
pylith::meshio::HDF5::isOpen(void) const
{ // isOpen
  return (_file == -1) ? false : true;
} // isOpen

// ----------------------------------------------------------------------
// Check if HDF5 file has group.
bool
pylith::meshio::HDF5::hasGroup(const char* name)
{ // hasGroup
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
    if (H5O_TYPE_GROUP == info.type)
      exists = true;
    err = H5Oclose(obj);
    assert(err >= 0);
  } // if

  PYLITH_METHOD_RETURN(exists);
} // hasGroup

// ----------------------------------------------------------------------
// Check if HDF5 file has dataset.
bool
pylith::meshio::HDF5::hasDataset(const char* name)
{ // hasDataset
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
    if (H5O_TYPE_DATASET == info.type)
      exists = true;
    err = H5Oclose(obj);
    assert(err >= 0);
  } // if

  PYLITH_METHOD_RETURN(exists);
} // hasDataset

// ----------------------------------------------------------------------
// Get topology metadata.
void
pylith::meshio::HDF5::getDatasetDims(hsize_t** dims,
				     int* ndims,
				     const char* parent,
				     const char* name)
{ // getDatasetDims
  PYLITH_METHOD_BEGIN;

  assert(dims);
  assert(ndims);
  assert(isOpen());

  try {
    // Open group
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t group = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
    hid_t group = H5Gopen(_file, parent);
#endif
    if (group < 0)
      throw std::runtime_error("Could not open group.");

    // Open the dataset
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
#else
    hid_t dataset = H5Dopen(group, name);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not open dataset.");

    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
      throw std::runtime_error("Could not get dataspace.");

    *ndims = H5Sget_simple_extent_ndims(dataspace);
    delete[] *dims; *dims = (*ndims > 0) ? new hsize_t[*ndims] : 0;
    H5Sget_simple_extent_dims(dataspace, *dims, 0);

    herr_t err = H5Sclose(dataspace);
    if (err < 0)
      throw std::runtime_error("Could not close dataspace.");

    err = H5Dclose(dataset);
    if (err < 0)
      throw std::runtime_error("Could not close dataset.");

    err = H5Gclose(group);
    if (err < 0)
      throw std::runtime_error("Could not close group.");

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

// ----------------------------------------------------------------------
// Get names of datasets in group.
void
pylith::meshio::HDF5::getGroupDatasets(string_vector* names,
				       const char* parent)
{ // getGroupDatasets
  PYLITH_METHOD_BEGIN;

  assert(names);
  assert(isOpen());

  try {
    // Open group
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t group = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
    hid_t group = H5Gopen(_file, parent);
#endif
    if (group < 0)
      throw std::runtime_error("Could not open group.");

    H5G_info_t ginfo;
    herr_t err = H5Gget_info(group, &ginfo);
    if (err < 0)
      throw std::runtime_error("Could not get group info.");
    const int gsize = ginfo.nlinks;

    names->resize(gsize);
    for (int i=0, index=0; i < gsize; ++i) {
      char buffer[256];
      ssize_t namelen = H5Lget_name_by_idx(group, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, buffer, 256, H5P_DEFAULT);assert(namelen > 0);
      (*names)[index++] = buffer;
    } // for

    err = H5Gclose(group);
    if (err < 0)
      throw std::runtime_error("Could not close group.");

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

// ----------------------------------------------------------------------
// Create group.
void
pylith::meshio::HDF5::createGroup(const char* name)
{ // createGroup
  PYLITH_METHOD_BEGIN;

  assert(name);

#if defined(PYLITH_HDF5_USE_API_18)
  hid_t group = H5Gcreate2(_file, name, 0, H5P_DEFAULT, H5P_DEFAULT);
#else // depracated HDF5 1.6 API
  hid_t group = H5Gcreate(_file, name, 0);
#endif
  if (group < 0) {
    std::ostringstream msg;
    msg << "Could not create group '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if

  herr_t err = H5Gclose(group);
  if (err < 0) {
    std::ostringstream msg;
    msg << "Could not close group '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if

  PYLITH_METHOD_END;
} // createGroup

// ----------------------------------------------------------------------
// Write scalar attribute.
void
pylith::meshio::HDF5::writeAttribute(const char* parent,
				     const char* name,
				     const void* value,
				     hid_t datatype)
{ // writeAttribute
  PYLITH_METHOD_BEGIN;

  HDF5::writeAttribute(_file, parent, name, value, datatype);

  PYLITH_METHOD_END;
} // writeAttribute

// ----------------------------------------------------------------------
// Write scalar attribute (external HDF5 handle).
void
pylith::meshio::HDF5::writeAttribute(hid_t h5,
				     const char* parent,
				     const char* name,
				     const void* value,
				     hid_t datatype)
{ // writeAttribute
  PYLITH_METHOD_BEGIN;

  assert(parent);
  assert(name);
  assert(value);

  try {
    hid_t dataspace = H5Screate(H5S_SCALAR);
    if (dataspace < 0)
      throw std::runtime_error("Could not create dataspace for");

#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen2(h5, parent, H5P_DEFAULT);
#else
    hid_t dataset = H5Dopen(h5, parent);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not open parent dataset for");

#if defined(PYLITH_HDF5_USE_API_18)
    hid_t attribute = H5Acreate2(dataset, name,
				datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT);
#else
    hid_t attribute = H5Acreate(dataset, name,
				datatype, dataspace, H5P_DEFAULT);
#endif
    if (attribute < 0)
      throw std::runtime_error("Could not create");

    hid_t err = H5Awrite(attribute, datatype, value);
    if (err < 0)
      throw std::runtime_error("Could not write");

    err = H5Aclose(attribute);
    if (err < 0)
      throw std::runtime_error("Could not close");

    err = H5Dclose(dataset);
    if (err < 0)
      throw std::runtime_error("Could not close dataset for");

    err = H5Sclose(dataspace);
    if (err < 0)
      throw std::runtime_error("Could not close dataspace for");

  } catch (std::exception& err) {
    std::ostringstream msg;
    msg << err.what() << " attribute '" << name << "' of '" << parent << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_END;
} // writeAttribute

// ----------------------------------------------------------------------
// Write string attribute.
void
pylith::meshio::HDF5::writeAttribute(const char* parent,
				     const char* name,
				     const char* value)
{ // writeAttribute
  PYLITH_METHOD_BEGIN;

  HDF5::writeAttribute(_file, parent, name, value);

  PYLITH_METHOD_END;
} // writeAttribute

// ----------------------------------------------------------------------
// Write string attribute (external handle to HDF5 file).
void
pylith::meshio::HDF5::writeAttribute(hid_t h5,
				     const char* parent,
				     const char* name,
				     const char* value)
{ // writeAttribute
  PYLITH_METHOD_BEGIN;

  assert(h5);
  assert(parent);
  assert(name);
  assert(value);

  try {
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen2(h5, parent, H5P_DEFAULT);
#else
    hid_t dataset = H5Dopen(h5, parent);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not open parent dataset for");

    hid_t dataspace = H5Screate(H5S_SCALAR);
    if (dataspace < 0)
      throw std::runtime_error("Could not create dataspace for");

    hid_t datatype = H5Tcopy(H5T_C_S1);
    if (datatype < 0)
      throw std::runtime_error("Could not create datatype for");

    herr_t err = H5Tset_size(datatype, strlen(value)+1);
    if (err < 0)
      throw std::runtime_error("Could not set size of");

#if defined(PYLITH_HDF5_USE_API_18)
    hid_t attribute = H5Acreate2(dataset, name,
				datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT);
#else
    hid_t attribute = H5Acreate(dataset, name,
				datatype, dataspace, H5P_DEFAULT);
#endif
    if (attribute < 0)
      throw std::runtime_error("Could not create");

    err = H5Awrite(attribute, datatype, value);
    if (err < 0)
      throw std::runtime_error("Could not write");

    err = H5Aclose(attribute);
    if (err < 0)
      throw std::runtime_error("Could not close");

    err = H5Tclose(datatype);
    if (err < 0)
      throw std::runtime_error("Could not close datatype for");

    err = H5Dclose(dataset);
    if (err < 0)
      throw std::runtime_error("Could not close dataset for");

    err = H5Sclose(dataspace);
    if (err < 0)
      throw std::runtime_error("Could not close dataspace for");

  } catch (std::exception& err) {
    std::ostringstream msg;
    msg << err.what() << " attribute '" << name << "' of '" << parent << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_END;
} // writeAttribute

// ----------------------------------------------------------------------
// Read scalar attribute.
void
pylith::meshio::HDF5::readAttribute(const char* parent,
				    const char* name,
				    void* value,
				    hid_t datatype)
{ // readAttribute
  PYLITH_METHOD_BEGIN;

  assert(parent);
  assert(name);
  assert(value);

  try {
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen2(_file, parent, H5P_DEFAULT);
#else
    hid_t dataset = H5Dopen(_file, parent);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not open parent dataset for");

    hid_t attribute = H5Aopen(dataset, name, H5P_DEFAULT);
    if (attribute < 0)
      throw std::runtime_error("Could not open");

    hid_t dtype = H5Aget_type(attribute);
    if (dtype < 0)
      throw std::runtime_error("Could not get datatype of");

    hid_t err = H5Aread(attribute, dtype, value);
    if (err < 0)
      throw std::runtime_error("Could not read");

    err = H5Tclose(dtype);
    if (err < 0)
      throw std::runtime_error("Could not close datatype for");

    err = H5Aclose(attribute);
    if (err < 0)
      throw std::runtime_error("Could not close");

    err = H5Dclose(dataset);
    if (err < 0)
      throw std::runtime_error("Could not close dataset for");

  } catch (std::exception& err) {
    std::ostringstream msg;
    msg << err.what() << " attribute '" << name << "' of '" << parent << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_END;
} // readAttribute

// ----------------------------------------------------------------------
// Read string attribute.
std::string
pylith::meshio::HDF5::readAttribute(const char* parent,
				    const char* name)
{ // readAttribute
  PYLITH_METHOD_BEGIN;

  assert(parent);
  assert(name);

  std::string value;

  try {
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen2(_file, parent, H5P_DEFAULT);
#else
    hid_t dataset = H5Dopen(_file, parent);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not open parent dataset for");

    hid_t attribute = H5Aopen(dataset, name, H5P_DEFAULT);
    if (attribute < 0)
      throw std::runtime_error("Could not open");

    hid_t datatype = H5Aget_type(attribute);
    if (datatype < 0)
      throw std::runtime_error("Could not get datatype of");

    // :TODO: Check that datatype is a string

    const int len = H5Tget_size(datatype);
    if (len <= 0)
      throw std::runtime_error("Nonpositive size for datatype of");

    char* buffer = (len > 0) ? new char[len] : 0;

    hid_t err = H5Aread(attribute, datatype, (void*)buffer);
    value = buffer;
    delete[] buffer; buffer = 0;

    if (err < 0)
      throw std::runtime_error("Could not read");

    err = H5Tclose(datatype);
    if (err < 0)
      throw std::runtime_error("Could not close datatype for");

    err = H5Aclose(attribute);
    if (err < 0)
      throw std::runtime_error("Could not close");

    err = H5Dclose(dataset);
    if (err < 0)
      throw std::runtime_error("Could not close dataset for");

  } catch (std::exception& err) {
    std::ostringstream msg;
    msg << err.what() << " attribute '" << name << "' of '" << parent << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_RETURN(std::string(value));
} // readAttribute

// ----------------------------------------------------------------------
// Create dataset.
void
pylith::meshio::HDF5::createDataset(const char* parent,
				    const char* name,
				    const hsize_t* maxDims,
				    const hsize_t* dimsChunk,
				    const int ndims,
				    hid_t datatype)
{ // createDataset
  PYLITH_METHOD_BEGIN;

  assert(parent);
  assert(name);
  assert(maxDims);
  assert(dimsChunk);

  try {
    // Open group
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t group = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
    hid_t group = H5Gopen(_file, parent);
#endif
    if (group < 0)
      throw std::runtime_error("Could not open group.");

    // Create the dataspace
    hsize_t *curDims = (ndims > 0) ? new hsize_t[ndims] : 0;
    if (ndims > 0)
      curDims[0] = (maxDims[0] == H5S_UNLIMITED) ? 1 : maxDims[0];
    for (int i=1; i < ndims; ++i)
      curDims[i] = maxDims[i];
    hid_t dataspace = H5Screate_simple(ndims, curDims, maxDims);
    delete[] curDims; curDims = 0;
    if (dataspace < 0)
      throw std::runtime_error("Could not create dataspace.");

    // Create chunked dataset
    hid_t property = H5Pcreate(H5P_DATASET_CREATE);
    if (property < 0)
      throw std::runtime_error("Could not create property for dataset.");

    herr_t err = H5Pset_chunk(property, ndims, dimsChunk);
    if (err < 0)
      throw std::runtime_error("Could not set chunk.");

    // Set gzip compression level for chunk.
    H5Pset_deflate(property, 6);

#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dcreate2(group, name,
			      datatype, dataspace, H5P_DEFAULT,
			      property, H5P_DEFAULT);
#else
    hid_t dataset = H5Dcreate(group, name,
			      datatype, dataspace, property);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not create dataset.");

    err = H5Dclose(dataset);
    if (err < 0)
      throw std::runtime_error("Could not close dataset.");

    err = H5Pclose(property);
    if (err < 0)
      throw std::runtime_error("Could not close property.");

    err = H5Sclose(dataspace);
    if (err < 0)
      throw std::runtime_error("Could not close dataspace.");

    err = H5Gclose(group);
    if (err < 0)
      throw std::runtime_error("Could not close group.");

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

// ----------------------------------------------------------------------
// Append slice to dataset.
void
pylith::meshio::HDF5::writeDatasetChunk(const char* parent,
					const char* name,
					const void* data,
					const hsize_t* dims,
					const hsize_t* dimsChunk,
					const int ndims,
					const int chunk,
					hid_t datatype)
{ // writeDatasetChunk
  PYLITH_METHOD_BEGIN;

  assert(parent);
  assert(name);
  assert(data);
  assert(dims);
  assert(_file > 0);

  try {
    // Select hyperslab in file
    hsize_t* count = (ndims > 0) ? new hsize_t[ndims] : 0;
    hsize_t* stride = (ndims > 0) ? new hsize_t[ndims] : 0;
    hsize_t* offset = (ndims > 0) ? new hsize_t[ndims] : 0;
    for (int i=0; i < ndims; ++i) {
      count[i] = 1;
      stride[i] = 1;
      offset[i] = 0;
    } // for
    offset[0] = chunk;

    // Open group
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t group = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
    hid_t group = H5Gopen(_file, parent);
#endif
    if (group < 0)
      throw std::runtime_error("Could not open group.");

    // Open the dataset
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
#else
    hid_t dataset = H5Dopen(group, name);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not open dataset.");

#if defined(PYLITH_HDF5_USE_API_18)
    herr_t err = H5Dset_extent(dataset, dims);
#else
    herr_t err = H5Dextend(dataset, dims);
#endif
    if (err < 0)
      throw std::runtime_error("Could not set dataset extent.");


    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
      throw std::runtime_error("Could not get dataspace.");

    hid_t chunkspace = H5Screate_simple(ndims, dimsChunk, 0);
    if (chunkspace < 0)
      throw std::runtime_error("Could not create chunk dataspace.");

    err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
			      offset, stride, count, dimsChunk);
    delete[] count; count = 0;
    delete[] stride; stride = 0;
    delete[] offset; offset = 0;
    if (err < 0)
      throw std::runtime_error("Could not select hyperslab.");

    err = H5Dwrite(dataset, datatype, chunkspace, dataspace,
		   H5P_DEFAULT, data);
    if (err < 0)
      throw std::runtime_error("Could not write data.");

    err = H5Sclose(chunkspace);
    if (err < 0)
      throw std::runtime_error("Could not close chunk dataspace.");

    err = H5Sclose(dataspace);
    if (err < 0)
      throw std::runtime_error("Could not close dataspace.");

    err = H5Dclose(dataset);
    if (err < 0)
      throw std::runtime_error("Could not close dataset.");

    err = H5Gclose(group);
    if (err < 0)
      throw std::runtime_error("Could not close group.");

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

// ----------------------------------------------------------------------
// Read dataset slice.
void
pylith::meshio::HDF5::readDatasetChunk(const char* parent,
				       const char* name,
				       char** const data,
				       hsize_t** const dimsChunk,
				       int* const ndims,
				       const int chunk,
				       hid_t datatype)
{ // readDatasetChunk
  PYLITH_METHOD_BEGIN;

  assert(parent);
  assert(name);
  assert(data);
  assert(dimsChunk);
  assert(_file > 0);

  try {
    // Open group
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t group = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
    hid_t group = H5Gopen(_file, parent);
#endif
    if (group < 0)
      throw std::runtime_error("Could not open group.");

    // Open the dataset
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
#else
    hid_t dataset = H5Dopen(group, name);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not open dataset.");

    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
      throw std::runtime_error("Could not get dataspace.");

    *ndims = H5Sget_simple_extent_ndims(dataspace);
    assert(*ndims > 0);
    hsize_t* dims = (*ndims > 0) ? new hsize_t[*ndims] : 0;
    H5Sget_simple_extent_dims(dataspace, dims, 0);

    // Select hyperslab in file
    delete[] *dimsChunk; *dimsChunk = (*ndims > 0) ? new hsize_t[*ndims] : 0;
    hsize_t* count = (*ndims > 0) ? new hsize_t[*ndims] : 0;
    hsize_t* stride = (*ndims > 0) ? new hsize_t[*ndims] : 0;
    hsize_t* offset = (*ndims > 0) ? new hsize_t[*ndims] : 0;

    for (int i=0; i < *ndims; ++i) {
      (*dimsChunk)[i] = dims[i];
      count[i] = 1;
      stride[i] = 1;
      offset[i] = 0;
    } // for
    (*dimsChunk)[0] = 1;
    offset[0] = chunk;

    hid_t chunkspace = H5Screate_simple(*ndims, *dimsChunk, 0);
    if (chunkspace < 0)
      throw std::runtime_error("Could not create chunk dataspace.");

    herr_t err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
				     offset, stride, count, *dimsChunk);
    delete[] count; count = 0;
    delete[] stride; stride = 0;
    delete[] offset; offset = 0;
    if (err < 0)
      throw std::runtime_error("Could not select hyperslab.");

    int sizeBytes = H5Tget_size(datatype);
    for (int i=0; i < *ndims; ++i)
      sizeBytes *= (*dimsChunk)[i];
    delete[] *data; *data = (sizeBytes > 0) ? new char[sizeBytes] : 0;
    delete[] dims; dims = 0;

    err = H5Dread(dataset, datatype, chunkspace, dataspace,
		  H5P_DEFAULT, (void*)*data);
    if (err < 0)
      throw std::runtime_error("Could not read data.");

    err = H5Sclose(chunkspace);
    if (err < 0)
      throw std::runtime_error("Could not close chunk dataspace.");

    err = H5Sclose(dataspace);
    if (err < 0)
      throw std::runtime_error("Could not close dataspace.");

    err = H5Dclose(dataset);
    if (err < 0)
      throw std::runtime_error("Could not close dataset.");

    err = H5Gclose(group);
    if (err < 0)
      throw std::runtime_error("Could not close group.");

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

// ----------------------------------------------------------------------
// Create dataset associated with data stored in a raw external binary
// file.
void
pylith::meshio::HDF5::createDatasetRawExternal(const char* parent,
					       const char* name,
					       const char* filename,
					       const hsize_t* maxDims,
					       const int ndims,
					       hid_t datatype)
{ // createDatasetRawExternal
  PYLITH_METHOD_BEGIN;

  assert(parent);
  assert(name);
  assert(filename);
  assert(maxDims);

  try {
    // Open group
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t group = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
    hid_t group = H5Gopen(_file, parent);
#endif
    if (group < 0)
      throw std::runtime_error("Could not open group.");

    // Create the dataspace
    hsize_t* curDims = (ndims > 0) ? new hsize_t[ndims] : 0;
    if (ndims > 0)
      curDims[0] = (maxDims[0] == H5S_UNLIMITED) ? 1 : maxDims[0];
    for (int i=1; i < ndims; ++i)
      curDims[i] = maxDims[i];
    hid_t dataspace = H5Screate_simple(ndims, curDims, maxDims);
    delete[] curDims; curDims = 0;
    if (dataspace < 0)
      throw std::runtime_error("Could not create dataspace.");

    // Create property for external dataset
    hid_t property = H5Pcreate(H5P_DATASET_CREATE);
    if (property < 0)
      throw std::runtime_error("Could not create property for dataset.");

    // Set external file
    const off_t offset = 0;
    herr_t err = H5Pset_external(property, filename, offset, H5F_UNLIMITED);
    if (err < 0)
      throw std::runtime_error("Could not set external file property.");

#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dcreate2(group, name,
			      datatype, dataspace, H5P_DEFAULT,
			      property, H5P_DEFAULT);
#else
    hid_t dataset = H5Dcreate(group, name,
			      datatype, dataspace, property);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not create dataset.");

    err = H5Dclose(dataset);
    if (err < 0)
      throw std::runtime_error("Could not close dataset.");

    err = H5Pclose(property);
    if (err < 0)
      throw std::runtime_error("Could not close property.");

    err = H5Sclose(dataspace);
    if (err < 0)
      throw std::runtime_error("Could not close dataspace.");

    err = H5Gclose(group);
    if (err < 0)
      throw std::runtime_error("Could not close group.");

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

// ----------------------------------------------------------------------
// Create dataset associated with data stored in a raw external binary
// file.
void
pylith::meshio::HDF5::extendDatasetRawExternal(const char* parent,
					       const char* name,
					       const hsize_t* dims,
					       const int ndims)
{ // extendDatasetRawExternal
  PYLITH_METHOD_BEGIN;

  assert(parent);
  assert(name);
  assert(dims);

  try {
    // Open group
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t group = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
    hid_t group = H5Gopen(_file, parent);
#endif
    if (group < 0)
      throw std::runtime_error("Could not open group.");

    // Open dataset.
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
#else
    hid_t dataset = H5Dopen(group, name);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not open dataset.");

#if defined(PYLITH_HDF5_USE_API_18)
    herr_t err = H5Dset_extent(dataset, dims);
#else
    herr_t err = H5Dextend(dataset, dims);
#endif
    if (err < 0)
      throw std::runtime_error("Could not set dataset extent.");

    err = H5Dclose(dataset);
    if (err < 0)
      throw std::runtime_error("Could not close dataset.");

    err = H5Gclose(group);
    if (err < 0)
      throw std::runtime_error("Could not close group.");

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

// ----------------------------------------------------------------------
// Write dataset comprised of an array of strings (external HDF5 handle).
void
pylith::meshio::HDF5::writeDataset(const char* parent,
				   const char* name,
				   const char* const* sarray,
				   const int nstrings)
{ // writeDataset
  PYLITH_METHOD_BEGIN;

  int slen = 0;
  for (int i=0; i < nstrings; ++i) {
      slen = std::max(slen, int(strlen(sarray[i])));
  } // for
  slen += 1; // add space for null terminator.

  char* strfixedlen = NULL;
  try {
      strfixedlen = (nstrings*slen > 0) ? new char[nstrings*slen] : 0;
      for (int i=0; i < nstrings; ++i) {
	  const int index = i*slen;
	  strncpy(&strfixedlen[index], sarray[i], slen-1);
	  strfixedlen[index+slen-1] = '\0';
	  // Fill rest of string with null.
	  for (int j=strlen(sarray[i]); j < slen; ++j) {
	      strfixedlen[index+j] = '\0';
	  } // for
      } // for

      HDF5::writeDataset(_file, parent, name, strfixedlen, nstrings, slen);

      delete[] strfixedlen; strfixedlen = NULL;
  } catch (const std::exception& err) {
      delete[] strfixedlen; strfixedlen = NULL;

      std::ostringstream msg;
      msg << "Error occurred while writing dataset '" << parent << "/" << name << "':\n" << err.what();
      throw std::runtime_error(msg.str());
  } catch (...) {
      delete[] strfixedlen; strfixedlen = NULL;

      std::ostringstream msg;
      msg << "Unknown error occurred while writing dataset '" << name << "'.";
      throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_END;
} // writeDataset

// ----------------------------------------------------------------------
// Write dataset comprised of an array of strings.
void
pylith::meshio::HDF5::writeDataset(const char* parent,
				   const char* name,
				   const char* sarray,
				   const int nstrings,
				   const int slen)
{ // writeDataset
  PYLITH_METHOD_BEGIN;

  HDF5::writeDataset(_file, parent, name, sarray, nstrings, slen);

  PYLITH_METHOD_END;
} // writeDataset

// ----------------------------------------------------------------------
// Write dataset comprised of an array of strings (external HDF5 handle).
void
pylith::meshio::HDF5::writeDataset(hid_t h5,
				   const char* parent,
				   const char* name,
				   const char* sarray,
				   const int nstrings,
				   const int slen)
{ // writeDataset
  PYLITH_METHOD_BEGIN;

  assert(h5);
  assert(parent);
  assert(name);
  assert(sarray);
  assert(nstrings > 0);
  assert(slen > 0);

  try {
    // Open group
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t group = H5Gopen2(h5, parent, H5P_DEFAULT);
#else
    hid_t group = H5Gopen(h5, parent);
#endif
    if (group < 0) throw std::runtime_error("Could not open group.");

    // Create the dataspace
    const int ndims = 1;
    hsize_t dims[ndims]; dims[0] = nstrings;
    hid_t dataspace = H5Screate_simple(ndims, dims, NULL);
    if (dataspace < 0) throw std::runtime_error("Could not create dataspace.");

    hid_t datatype = H5Tcopy(H5T_C_S1);
    if (datatype < 0) throw std::runtime_error("Could not create datatype.");
    herr_t err = H5Tset_size(datatype, slen);
    if (err < 0) throw std::runtime_error("Could not set size of datatype.");

#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dcreate2(group, name, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
    hid_t dataset = H5Dcreate(group, name, datatype, dataspace, H5P_DEFAULT);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not create dataset.");

    err = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, sarray);
    if (err < 0) throw std::runtime_error("Could not write dataset.");

    err = H5Dclose(dataset);
    if (err < 0) throw std::runtime_error("Could not close dataset.");

    err = H5Tclose(datatype);
    if (err < 0) throw std::runtime_error("Could not close datatype.");

    err = H5Sclose(dataspace);
    if (err < 0) throw std::runtime_error("Could not close dataspace.");

    err = H5Gclose(group);
    if (err < 0) throw std::runtime_error("Could not close group.");

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


// ----------------------------------------------------------------------
// Read dataset comprised of an array of strings.
pylith::string_vector
pylith::meshio::HDF5::readDataset(const char* parent,
				  const char* name)
{ // readDataset
  PYLITH_METHOD_BEGIN;

  assert(_file);
  assert(parent);
  assert(name);

  pylith::string_vector data;
  char* strfixedlen = 0;
  try {
    // Open group
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t group = H5Gopen2(_file, parent, H5P_DEFAULT);
#else
    hid_t group = H5Gopen(_file, parent);
#endif
    if (group < 0)
      throw std::runtime_error("Could not open group.");

    // Open the dataset
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
#else
    hid_t dataset = H5Dopen(group, name);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not open dataset.");

    // Get the datatype
    hid_t datatype = H5Dget_type(dataset);
    if (datatype < 0)
      throw std::runtime_error("Could not get datatype.");
    const int slen = H5Tget_size(datatype);
    if (slen < 0)
      throw std::runtime_error("Could not get size of datatype.");

    // Get the dataspace
    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
      throw std::runtime_error("Could not get dataspace.");
    const int ndims = H5Sget_simple_extent_ndims(dataspace);
    if (ndims != 1)
      throw std::runtime_error("Expected 1 dimension for string dataset.");
    hsize_t dims[1];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    const int nstrings = dims[0];
    if (nstrings <= 0)
      throw std::runtime_error("Zero size for dataset.");

    strfixedlen = (nstrings*slen > 0) ? new char[nstrings*slen] : 0;

    // Read the dataset
    herr_t err = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, strfixedlen);

    data.resize(nstrings);
    for (int i=0; i < nstrings; ++i) {
      data[i] = &strfixedlen[i*slen];
    } // for
    delete[] strfixedlen; strfixedlen = 0;

    err = H5Dclose(dataset);
    if (err < 0)
      throw std::runtime_error("Could not close dataset.");

    err = H5Tclose(datatype);
    if (err < 0)
      throw std::runtime_error("Could not close datatype.");

    err = H5Sclose(dataspace);
    if (err < 0)
      throw std::runtime_error("Could not close dataspace.");

    err = H5Gclose(group);
    if (err < 0)
      throw std::runtime_error("Could not close group.");

  } catch (const std::exception& err) {
    delete[] strfixedlen; strfixedlen = 0;

    std::ostringstream msg;
    msg << "Error occurred while creating dataset '"
	<< parent << "/" << name << "':\n"
	<< err.what();
    throw std::runtime_error(msg.str());
  } catch (...) {
    delete[] strfixedlen; strfixedlen = 0;

    std::ostringstream msg;
    msg << "Unknown error occurred while writing dataset '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_RETURN(data);
} // readDataset


// End of file
