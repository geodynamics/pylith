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

#include "ExodusII.hh" // implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "petsc.h" // USES MPI_Comm

// :KLUDGE: Prevent NetCDF from definining MPI types
#define MPI_INCLUDED
#include <netcdf.h> // USES netcdf

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::ExodusII::ExodusII(void) :
  _filename(""),
  _file(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Constructor with filename
pylith::meshio::ExodusII::ExodusII(const char* filename) :
  _filename(filename),
  _file(0)
{ // constructor
  open();
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::ExodusII::~ExodusII(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate data structures.
void
pylith::meshio::ExodusII::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  close();

  PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set filename.
void
pylith::meshio::ExodusII::filename(const char* name)
{ // filename
  _filename = name;
} // filename

// ----------------------------------------------------------------------
// Get filename.
const char*
pylith::meshio::ExodusII::filename(void) const
{ // filename
  return _filename.c_str();
} // filename

// ----------------------------------------------------------------------
// Open file.
void
pylith::meshio::ExodusII::open(void)
{ // open
  PYLITH_METHOD_BEGIN;

  close();

  const int err = nc_open(_filename.c_str(), NC_NOWRITE, &_file);
  if (err != NC_NOERR) {
    std::ostringstream msg;
    msg << "Could not open ExodusII file '" << _filename << "' for reading.\n";
    throw std::runtime_error(msg.str());
  } // if

  PYLITH_METHOD_END;
} // open

// ----------------------------------------------------------------------
// Close file.
void
pylith::meshio::ExodusII::close(void)
{ // close
  PYLITH_METHOD_BEGIN;

  if (_file) {
    int err = nc_close(_file);_file=0;
    if (err != NC_NOERR) {
      std::ostringstream msg;
      msg << "Could not close ExodusII file '" << _filename << "'.\n";
      throw std::runtime_error(msg.str());
    } // if
  } // if

  PYLITH_METHOD_END;
} // close

// ----------------------------------------------------------------------
// Check if Cubit Exodus file constains dimension.
bool
pylith::meshio::ExodusII::hasDim(const char* name,
				 int* id) const
{ // hasDim
  PYLITH_METHOD_BEGIN;

  assert(_file);

  bool found = false;

  int i = -1;
  int err = nc_inq_dimid(_file, name, &i);
  if (err == NC_NOERR) {
    found = true;
    if (id) {
      *id = i;
    } // if
  } else if (id) {
    *id = -1;
  } // if/else
  
  PYLITH_METHOD_RETURN(found);
} // hasDim

// ----------------------------------------------------------------------
// Check if Cubit Exodus file constains attribute.
bool
pylith::meshio::ExodusII::hasAtt(const char* name,
				 int* id) const
{ // hasAtt
  PYLITH_METHOD_BEGIN;

  assert(_file);

  bool found = false;

  int i = -1;
  int err = nc_inq_attid(_file, NC_GLOBAL, name, &i);
  if (err == NC_NOERR) {
    found = true;
    if (id) {
      *id = i;
    } // if
  } else if (id) {
    *id = -1;
  } // if/else
  
  PYLITH_METHOD_RETURN(found);
} // hasAtt

// ----------------------------------------------------------------------
// Check if Cubit Exodus file constains variable.
bool
pylith::meshio::ExodusII::hasVar(const char* name,
				 int* id) const
{ // hasVar
  PYLITH_METHOD_BEGIN;

  assert(_file);

  bool found = false;

  int i = -1;
  int err = nc_inq_varid(_file, name, &i);
  if (err == NC_NOERR) {
    found = true;
    if (id) {
      *id = i;
    } // if
  } else if(id) {
    *id = -1;
  } // if/else

  PYLITH_METHOD_RETURN(found);
} // hasVar

// ----------------------------------------------------------------------
// Get value for dimension.
int
pylith::meshio::ExodusII::getDim(const char* name) const
{ // getDim
  PYLITH_METHOD_BEGIN;

  assert(_file);

  int id = -1;
  if (!hasDim(name, &id)) {
    std::ostringstream msg;
    msg << "Missing dimension '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if

  size_t dimSize = 0;
  int err = nc_inq_dimlen(_file, id, &dimSize);
  if (err != NC_NOERR) {
    std::ostringstream msg;
    msg << "Could not get dimension '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  PYLITH_METHOD_RETURN(dimSize);
} // getDim

// ----------------------------------------------------------------------
// Get values for variable as an array of PylithScalars.
void
pylith::meshio::ExodusII::getVar(PylithScalar* values,
				 int* dims,
				 int ndims,
				 const char* name) const
{ // getVar
  PYLITH_METHOD_BEGIN;

  assert(_file);
  assert(values);

  int vid = -1;
  if (!hasVar(name, &vid)) {
    std::ostringstream msg;
    msg << "Missing real variable '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if

  int vndims = 0;
  int err = nc_inq_varndims(_file, vid, &vndims);
  if (ndims != vndims) {
    std::ostringstream msg;
    msg << "Expecting " << ndims << " dimensions for variable '" << name
	<< "' but variable only has " << vndims << " dimensions.";
    throw std::runtime_error(msg.str());
  } // if

  int* dimIds = (ndims > 0) ? new int[ndims] : 0;
  err = nc_inq_vardimid(_file, vid, dimIds);
  if (err != NC_NOERR) {
    delete[] dimIds; dimIds = 0;
    std::ostringstream msg;
    msg << "Could not get dimensions for variable '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  for (int iDim=0; iDim < ndims; ++iDim) {
    size_t dimSize = 0;
    err = nc_inq_dimlen(_file, dimIds[iDim], &dimSize);
    if (err != NC_NOERR) {
      delete[] dimIds; dimIds = 0;
      std::ostringstream msg;
      msg << "Could not get dimension '" << iDim << "' for variable '" << name << "'.";
      throw std::runtime_error(msg.str());
    } // if
    if (size_t(dims[iDim]) != dimSize) {
      delete[] dimIds; dimIds = 0;
      std::ostringstream msg;
      msg << "Expecting dimension " << iDim << " of variable '" << name
	  << "' to be " << dims[iDim] << ", but dimension is " << dimSize
	  << ".";
      throw std::runtime_error(msg.str());
    } // if
  } // for
  delete[] dimIds; dimIds = 0;

  if (sizeof(PylithScalar) == sizeof(double)) {
    err = nc_get_var_double(_file, vid, values);
  } else {
    assert(0);
    throw std::logic_error("Unknown size of PylithScalar in ExodusII::getVar().");
  } // if/else
  if (err != NC_NOERR) {
    std::ostringstream msg;
    msg << "Coult not get values for variable '" << name << ".";
    throw std::runtime_error(msg.str());
  } // if

  PYLITH_METHOD_END;
} // getVar

// ----------------------------------------------------------------------
// Get values for variable as an array of ints.
void
pylith::meshio::ExodusII::getVar(int* values,
				 int* dims,
				 int ndims,
				 const char* name) const
{ // getVar
  PYLITH_METHOD_BEGIN;

  assert(_file);
  assert(values);

  int vid = -1;
  if (!hasVar(name, &vid)) {
    std::ostringstream msg;
    msg << "Missing integer variable '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if

  int vndims = 0;
  int err = nc_inq_varndims(_file, vid, &vndims);
  if (ndims != vndims) {
    std::ostringstream msg;
    msg << "Expecting " << ndims << " dimensions for variable '" << name
	<< "' but variable only has " << vndims << " dimensions.";
    throw std::runtime_error(msg.str());
  } // if

  int* dimIds = (ndims > 0) ? new int[ndims] : 0;
  err = nc_inq_vardimid(_file, vid, dimIds);
  if (err != NC_NOERR) {
    delete[] dimIds; dimIds = 0;
    std::ostringstream msg;
    msg << "Could not get dimensions for variable '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  for (int iDim=0; iDim < ndims; ++iDim) {
    size_t dimSize = 0;
    err = nc_inq_dimlen(_file, dimIds[iDim], &dimSize);
    if (err != NC_NOERR) {
      delete[] dimIds; dimIds = 0;
      std::ostringstream msg;
      msg << "Could not get dimension '" << iDim << "' for variable '" << name << "'.";
      throw std::runtime_error(msg.str());
    } // if
    if (size_t(dims[iDim]) != dimSize) {
      delete[] dimIds; dimIds = 0;
      std::ostringstream msg;
      msg << "Expecting dimension " << iDim << " of variable '" << name
	  << "' to be " << dims[iDim] << ", but dimension is " << dimSize
	  << ".";
      throw std::runtime_error(msg.str());
    } // if
  } // for
  delete[] dimIds; dimIds = 0;

  err = nc_get_var_int(_file, vid, values);
  if (err != NC_NOERR) {
    std::ostringstream msg;
    msg << "Coult not get values for variable '" << name << ".";
    throw std::runtime_error(msg.str());
  } // if

  PYLITH_METHOD_END;
} // getVar

// ----------------------------------------------------------------------
// Get values for variable as an array of strings.
void
pylith::meshio::ExodusII::getVar(string_vector* values,
				 int dim,
				 const char* name) const
{ // getVar
  PYLITH_METHOD_BEGIN;

  assert(_file);
  assert(values);
  
  int vid = -1;
  if (!hasVar(name, &vid)) {
    std::ostringstream msg;
    msg << "Missing string variable '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if

  int vndims = 0;
  const int ndims = 2;
  int err = nc_inq_varndims(_file, vid, &vndims);
  if (ndims != vndims) {
    std::ostringstream msg;
    msg << "Expecting " << ndims << " dimensions for variable '" << name
	<< "' but variable only has " << vndims << " dimensions.";
    throw std::runtime_error(msg.str());
  } // if

  int dimIds[ndims];
  err = nc_inq_vardimid(_file, vid, dimIds);
  if (err != NC_NOERR) {
    std::ostringstream msg;
    msg << "Could not get dimensions for variable '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if

  int iDim = 0;
  size_t dimSize = 0;
  err = nc_inq_dimlen(_file, dimIds[0], &dimSize);
  if (err != NC_NOERR) {
    std::ostringstream msg;
    msg << "Could not get dimension '" << iDim << "' for variable '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if
  if (size_t(dim) != dimSize) {
    std::ostringstream msg;
    msg << "Expecting dimension " << iDim << " of variable '" << name
	<< "' to be " << dim << ", but dimension is " << dimSize
	<< ".";
    throw std::runtime_error(msg.str());
  } // if

  iDim = 1;
  err = nc_inq_dimlen(_file, dimIds[1], &dimSize);
  if (err != NC_NOERR) {
    std::ostringstream msg;
    msg << "Could not get dimension '" << iDim << "' for variable '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if
  const size_t bufferSize = dimSize;
  char* buffer = (bufferSize > 0) ? new char[bufferSize] : 0;
  size_t indices[2] = { 0, 0 };
  size_t chunk[2] = { 1, bufferSize };
  for (int i=0; i < dim; ++i) {
    indices[0] = i;
    err = nc_get_vara_text(_file, vid, indices, chunk, buffer);
    if (err != NC_NOERR) {
      delete[] buffer; buffer = 0;
      std::ostringstream msg;
      msg << "Could not read string in variable '" << name << "'.";
      throw std::runtime_error(msg.str());
    } // if
    (*values)[i] = buffer;
    //std::cout << "GROUP: '" << (*values)[i] << "'." << std::endl;
  } // for
  delete[] buffer; buffer = 0;

  PYLITH_METHOD_END;
} // getVar


// End of file 
