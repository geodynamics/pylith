// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include "HDF5.hh" // implementation of class methods

#include <cstring> // USES strlen()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 8 && !defined(H5_USE_16_API)
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
			   hid_t mode,
			   const bool create)
{ // constructor
  if (create) {
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
			   hid_t mode,
			   const bool create)
{ // open
  if (_file >= 0) {
    throw std::runtime_error("HDF5 file already open.");
  } // if

  if (create) {
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
} // constructor

// ----------------------------------------------------------------------
// Close HDF5 file.
void
pylith::meshio::HDF5::close(void)
{ // close
  if (_file >= 0) {
    herr_t err = H5Fclose(_file);
    if (err < 0) 
      throw std::runtime_error("Could not close HDF5 file.");
  } // if
  _file = -1;
} // close

// ----------------------------------------------------------------------
// Create group.
void
pylith::meshio::HDF5::createGroup(const char* name)
{ // createGroup
#if defined(PYLITH_HDF5_USE_API_18)
  hid_t group = H5Gcreate(_file, name, 0, H5P_DEFAULT, H5P_DEFAULT);
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
} // createGroup

// ----------------------------------------------------------------------
// Write scalar attribute.
void
pylith::meshio::HDF5::writeAttribute(const char* parent,
				     const char* name,
				     const void* value,
				     hid_t datatype)
{ // writeAttribute
  try {
    hid_t dataspace = H5Screate(H5S_SCALAR);
    if (dataspace < 0)
      throw std::runtime_error("Could not create dataspace for");

#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen(_file, parent, H5P_DEFAULT);
#else
    hid_t dataset = H5Dopen(_file, parent);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not open parent dataset for");

#if defined(PYLITH_HDF5_USE_API_18)
    hid_t attribute = H5Acreate(dataset, name,
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
} // writeAttribute

// ----------------------------------------------------------------------
// Read scalar attribute.
void
pylith::meshio::HDF5::readAttribute(const char* parent,
				    const char* name,
				    void* value,
				    hid_t datatype)
{ // readAttribute
  try {
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen(_file, parent, H5P_DEFAULT);
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

    if (H5Tequal(dtype, datatype) <= 0)
      throw std::runtime_error("Wrong datatype specified for");

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
} // readAttribute

// ----------------------------------------------------------------------
// Write string attribute.
void
pylith::meshio::HDF5::writeAttribute(const char* parent,
				     const char* name,
				     const char* value)
{ // writeAttribute
  try {
    hid_t dataspace = H5Screate(H5S_SCALAR);
    if (dataspace < 0) 
      throw std::runtime_error("Could not create dataspace for");

#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen(_file, parent, H5P_DEFAULT);
#else
    hid_t dataset = H5Dopen(_file, parent);
#endif
    if (dataset < 0) 
      throw std::runtime_error("Could not open parent dataset for");

    hid_t datatype = H5Tcopy(H5T_C_S1);
    if (datatype < 0) 
      throw std::runtime_error("Could not create datatype for");

    herr_t err = H5Tset_size(datatype, strlen(value)+1);
    if (err < 0) 
      throw std::runtime_error("Could not set size of");

#if defined(PYLITH_HDF5_USE_API_18)
    hid_t attribute = H5Acreate(dataset, name,
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
} // writeAttribute

// ----------------------------------------------------------------------
// Read string attribute.
const char*
pylith::meshio::HDF5::readAttribute(const char* parent,
				    const char* name)
{ // readAttribute
  std::string value;

  try {
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen(_file, parent, H5P_DEFAULT);
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

  return value.c_str();
} // readAttribute

// ----------------------------------------------------------------------
// Create dataset.
void
pylith::meshio::HDF5::createDataset(const char* parent,
				    const char* name,
				    const hsize_t* dims,
				    const hsize_t ndims,
				    hid_t datatype)
{ // createDataset
  try {
    // Open group
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t group = H5Gopen(_file, parent, H5P_DEFAULT);
#else
    hid_t group = H5Gopen(_file, parent);
#endif
    if (group < 0) 
      throw std::runtime_error("Could not open group.");

    // Create the dataspace
    hid_t dataspace = H5Screate_simple(ndims, dims, 0);
    if (dataspace < 0)
      throw std::runtime_error("Could not create dataspace.");
      
    // Create chunked dataset
    hid_t property = H5Pcreate(H5P_DATASET_CREATE);
    if (property < 0)
      throw std::runtime_error("Could not create property for dataset.");

    herr_t err = H5Pset_chunk(property, ndims, dims);
    if (err < 0)
      throw std::runtime_error("Could not set chunk.");

    // Set gzip compression level for chunk.
    //H5Pset_deflate(property, 6);

#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dcreate(group, name,
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
    msg << "Unknown  occurred while creating dataset '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch
} // createDataset

// ----------------------------------------------------------------------
// Create dataset associated with data stored in a raw external binary
// file.
void
pylith::meshio::HDF5::createDatasetRawExternal(const char* parent,
					       const char* name,
					       const char* filename,
					       const hsize_t* dims,
					       const hsize_t ndims,
					       hid_t datatype)
{ // createDatasetRawExternal
  try {
    // Open group
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t group = H5Gopen(_file, parent, H5P_DEFAULT);
#else
    hid_t group = H5Gopen(_file, parent);
#endif
    if (group < 0) 
      throw std::runtime_error("Could not open group.");

    // Create the dataspace
    hid_t dataspace = H5Screate_simple(ndims, dims, 0);
    if (dataspace < 0)
      throw std::runtime_error("Could not create dataspace.");
      
    // Create chunked dataset
    hid_t property = H5Pcreate(H5P_DATASET_CREATE);
    if (property < 0)
      throw std::runtime_error("Could not create property for dataset.");

#if 0 // Not allowed for external storage?
    herr_t err = H5Pset_chunk(property, ndims, dims);
    if (err < 0)
      throw std::runtime_error("Could not set chunk.");

    // Set gzip compression level for chunk.
    //H5Pset_deflate(property, 6);
#endif

    // Set external file
    const off_t offset = 0;
    hsize_t sizeBytes = 0;
    for (int i=0; i < ndims; ++i)
      sizeBytes += dims[i] * H5Tget_size(datatype);
    herr_t err = H5Pset_external(property, filename, offset, sizeBytes);
    if (err < 0)
      throw std::runtime_error("Could not set external file property.");

#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dcreate(group, name,
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
    msg << "Unknown  occurred while creating dataset '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch
} // createDatasetRawExternal

// ----------------------------------------------------------------------
// Append slice to dataset.
void
pylith::meshio::HDF5::writeDatasetSlice(const char* parent,
					const char* name,
					const void* data,
					const hsize_t* dims,
					const hsize_t ndims,
					const int islice,
					hid_t datatype)
{ // writeDatasetSlice
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
    offset[0] = islice;

    // Open group
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t group = H5Gopen(_file, parent, H5P_DEFAULT);
#else
    hid_t group = H5Gopen(_file, parent);
#endif
    if (group < 0)
      throw std::runtime_error("Could not open group.");
    
    // Open the dataset
#if defined(PYLITH_HDF5_USE_API_18)
    hid_t dataset = H5Dopen(group, name, H5P_DEFAULT);
#else
    hid_t dataset = H5Dopen(group, name);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not open dataset.");
    
    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
      throw std::runtime_error("Could not get dataspace.");

    hid_t chunkspace = H5Screate_simple(ndims, dims, 0);
    if (chunkspace < 0)
      throw std::runtime_error("Could not create chunk dataspace.");

    herr_t err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
				     offset, stride, count, dims);
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
    msg << "Unknown  occurred while writing dataset '"
	<< parent << "/" << name << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch
} // writeDatasetSlice


// End of file
