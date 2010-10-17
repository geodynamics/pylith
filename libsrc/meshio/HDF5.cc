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

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::HDF5::HDF5(const char* filename,
			   hid_t mode)
{ // constructor
  _file = H5Fopen(filename, mode, H5P_DEFAULT);
  if (_file < 0) {
    std::ostringstream msg;
    msg << "Could not open HDF5 file '" << filename << "'.";
    throw std::runtime_error(msg.str());
  } // if
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::HDF5::~HDF5(void)
{ // destructor
  H5Fclose(_file);
} // destructor

// ----------------------------------------------------------------------
// Create group.
hid_t
pylith::meshio::HDF5::createGroup(const char* name)
{ // createGroup
  hid_t group = H5Gcreate(_file, name, 0, H5P_DEFAULT, H5P_DEFAULT);
  if (group < 0) {
    std::ostringstream msg;
    msg << "Coule not create group '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if

  return group;
} // createGroup

// ----------------------------------------------------------------------
// Write scalar attribute.
void
pylith::meshio::HDF5::writeAttribute(hid_t parent,
				     const char* name,
				     const void* value,
				     hid_t datatype)
{ // writeAttribute
  try {
    hid_t dataspace = H5Screate(H5S_SCALAR);
    if (dataspace < 0) {
      throw std::runtime_error("Could not create dataspace for");
    } // if
    hid_t attribute = H5Acreate(parent, name,
				datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    if (attribute < 0) {
      throw std::runtime_error("Could not create");
    } // if
    hid_t err = H5Awrite(attribute, datatype, value);
    if (err < 0) {
      throw std::runtime_error("Could not write");
    } // if
    err = H5Aclose(attribute);
    if (err < 0) {
      throw std::runtime_error("Could not close");
    } // if
  } catch (std::exception& err) {
    std::ostringstream msg;
    msg << err.what() << " attribute '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch
} // writeAttribute

// ----------------------------------------------------------------------
// Write string attribute.
void
pylith::meshio::HDF5::writeAttribute(hid_t parent,
				     const char* name,
				     const char* value)
{ // writeAttribute
  try {
    hid_t dataspace = H5Screate(H5S_SCALAR);
    if (dataspace < 0) {
      throw std::runtime_error("Could not create dataspace for");
    } // if
    hid_t attribute = H5Acreate(parent, name,
				H5T_C_S1, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    if (attribute < 0) {
      throw std::runtime_error("Could not create");
    } // if
    hid_t err = H5Awrite(attribute, H5T_C_S1, value);
    if (err < 0) {
      throw std::runtime_error("Could not write");
    } // if
    err = H5Aclose(attribute);
    if (err < 0) {
      throw std::runtime_error("Could not close");
    } // if
  } catch (std::exception& err) {
    std::ostringstream msg;
    msg << err.what() << " attribute '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch
} // writeAttribute

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
    // Create the dataspace
    hid_t dataspace = H5Screate_simple(ndims, dims, 0);
    if (dataspace < 0)
      throw std::runtime_error("Could not create dataspace.");
      
    // Create chunked dataset
    hid_t property = H5Pcreate(H5P_DATASET_CREATE);
    if (property < 0)
      throw std::runtime_error("Could not create property for dataset.");

    H5Pset_chunk(property, ndims, dims);
    //H5Pset_deflate(property, 6);

    std::string fullname = 
      std::string(parent) + std::string("/") + std::string(name);
#if defined(OLD_H5INTERFACE)
    hid_t dataset = H5Dcreate(_file, fullname.c_str(),
			      datatype, dataspace, property);
#else
    hid_t dataset = H5Dcreate(_file, fullname.c_str(),
			      datatype, dataspace, H5P_DEFAULT,
			      property, H5P_DEFAULT);
#endif
    if (dataset < 0)
      throw std::runtime_error("Could not create dataset.");
    H5Pclose(property);
    H5Sclose(dataspace);
    H5Dclose(dataset);
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error occurred while creating dataset '"
	<< parent << "/" << name << "':\n"
	<< err.what();
    throw std::runtime_error(msg.str());
  } catch (...) {
    std::ostringstream msg;
    msg << "Unknown  occurred while creating dataset '"
	<< parent << "/" << name << "'.";
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
    hid_t group = H5Gopen(_file, parent, H5P_DEFAULT);
    if (group < 0)
      throw std::runtime_error("Could not open group.");
    
    // Open the dataset
    hid_t dataset = H5Dopen(group, name, H5P_DEFAULT);
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
