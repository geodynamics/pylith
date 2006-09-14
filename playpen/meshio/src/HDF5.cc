// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include "HDF5.hh" // implementation of class methods

extern "C" {
#include "hdf5.h" // USES hdf5
}

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::HDF5::HDF5(const char* filename,
			   hid_t mode)
{ // constructor
  _file = H5Fopen(filename, mode, H5P_DEFAULT);
  if (_file < 0) {
    std::ostringstream msg;
    msg << "Could not open HDF5 mesh file '" << filename << "'.";
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
  hid_t group = H5Gcreate(_file, name, 0);
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
				     const void* pValue,
				     hid_t datatype)
{ // writeAttribute
  try {
    hid_t dataspace = H5Screate(H5S_SCALAR);
    if (dataspace < 0)
      throw std::runtime_error("Could not create dataspace for");
    hid_t attribute = H5Acreate(parent, name,
				datatype, dataspace, H5P_DEFAULT);
    if (attribute < 0)
      throw std::runtime_error("Could not create");
    hid_t err = H5Awrite(attribute, datatype, pValue);
    if (err < 0)
      throw std::runtime_error("Could not write");
    err = H5Aclose(attribute);
    if (err < 0)
      throw std::runtime_error("Could not close");
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
} // writeAttribute

// End of file 
