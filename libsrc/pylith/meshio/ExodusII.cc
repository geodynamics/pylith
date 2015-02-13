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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "ExodusII.hh" // implementation of class methods

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include "petsc.h" // USES MPI_Comm

// :KLUDGE: Prevent NetCDF from definining MPI types
#define MPI_INCLUDED
#include <netcdfcpp.h> // USES netcdf

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

  _file = new NcFile(_filename.c_str());
  if (!_file->is_valid()) {
    std::ostringstream msg;
    msg << "Could not open ExodusII file '" << _filename
	<< "' for reading.\n";
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

  if (_file)
    _file->close();
  delete _file; _file = 0;

  PYLITH_METHOD_END;
} // close

// ----------------------------------------------------------------------
// Check if Cubit Exodus file constains dimension.
bool
pylith::meshio::ExodusII::hasDim(const char* name) const
{ // hasDim
  PYLITH_METHOD_BEGIN;

  assert(_file);

  bool found = false;
  
  const int ndims = _file->num_dims();
  for (int i=0; i < ndims; ++i) {
    NcDim* dim = _file->get_dim(i);
    assert(dim);
    if (0 == strcmp(dim->name(), name)) {
      found = true;
      break;
    } // if
  } // for
  
  PYLITH_METHOD_RETURN(found);
} // hasDim

// ----------------------------------------------------------------------
// Check if Cubit Exodus file constains attribute.
bool
pylith::meshio::ExodusII::hasAtt(const char* name) const
{ // hasAtt
  PYLITH_METHOD_BEGIN;

  assert(_file);

  bool found = false;

  const int natts = _file->num_atts();
  for (int i=0; i < natts; ++i) {
    NcAtt* att = _file->get_att(i);
    assert(att);
    if (0 == strcmp(att->name(), name)) {
      found = true;
      delete att; att = 0;
      break;
    } // if
    delete att; att = 0;
  } // for
  
  PYLITH_METHOD_RETURN(found);
} // hasAtt

// ----------------------------------------------------------------------
// Check if Cubit Exodus file constains variable.
bool
pylith::meshio::ExodusII::hasVar(const char* name) const
{ // hasVar
  PYLITH_METHOD_BEGIN;

  assert(_file);

  bool found = false;

  const int nvars = _file->num_vars();
  for (int i=0; i < nvars; ++i) {
    NcVar* var = _file->get_var(i);
    assert(var);
    if (0 == strcmp(var->name(), name)) {
      found = true;
      break;
    } // if
  } // for

  PYLITH_METHOD_RETURN(found);
} // hasVar

// ----------------------------------------------------------------------
// Get value for dimension.
int
pylith::meshio::ExodusII::getDim(const char* name) const
{ // getDim
  PYLITH_METHOD_BEGIN;

  assert(_file);

  if (!hasDim(name)) {
    std::ostringstream msg;
    msg << "Missing dimension '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  NcDim* dim = _file->get_dim(name);
  assert(dim);
  PYLITH_METHOD_RETURN(dim->size());
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

  if (!hasVar(name)) {
    std::ostringstream msg;
    msg << "Missing real variable '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if

  NcVar* var = _file->get_var(name);
  assert(var);
  if (ndims != var->num_dims()) {
    std::ostringstream msg;
    msg << "Expecting " << ndims << " dimensions for variable '" << name
	<< "' but variable only has " << var->num_dims() << " dimensions.";
    throw std::runtime_error(msg.str());
  } // if
  
  for (int iDim=0; iDim < ndims; ++iDim) {
    NcDim* d = var->get_dim(iDim);
    assert(d);
    if (dims[iDim] != d->size()) {
      std::ostringstream msg;
      msg << "Expecting dimension " << iDim << " of variable '" << name
	  << "' to be " << dims[iDim] << ", but dimension is " << d->size()
	  << ".";
      throw std::runtime_error(msg.str());
    } // if
  } // for

  long* counts = var->edges();
  bool ok = var->get(values, counts);
  delete[] counts; counts = 0;
  if (!ok) {
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

  if (!hasVar(name)) {
    std::ostringstream msg;
    msg << "Missing integer variable '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if

  NcVar* var = _file->get_var(name);
  assert(var);
  if (ndims != var->num_dims()) {
    std::ostringstream msg;
    msg << "Expecting " << ndims << " dimensions for variable '" << name
	<< "' but variable only has " << var->num_dims() << " dimensions.";
    throw std::runtime_error(msg.str());
  } // if
  
  for (int iDim=0; iDim < ndims; ++iDim) {
    NcDim* d = var->get_dim(iDim);
    assert(d);
    if (dims[iDim] != d->size()) {
      std::ostringstream msg;
      msg << "Expecting dimension " << iDim << " of variable '" << name
	  << "' to be " << dims[iDim]<< ", but dimension is " << d->size()
	  << ".";
      throw std::runtime_error(msg.str());
    } // if
  } // for

  long* counts = var->edges();
  bool ok = var->get(values, counts);
  delete[] counts; counts = 0;
  if (!ok) {
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
  
  if (!hasVar(name)) {
    std::ostringstream msg;
    msg << "Missing string variable '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if

  NcVar* var = _file->get_var(name);
  assert(var);

  if (var->num_dims() != 2) {
    std::ostringstream msg;
    msg << "Expecting 2 dimensions for variable '" << name
	<< "' but variable only has " << var->num_dims() << " dimensions.";
    throw std::runtime_error(msg.str());
  } // if

  long* counts = var->edges();
  const int bufferSize = counts[1];
  char* buffer = (bufferSize > 0) ? new char[bufferSize] : 0;
  for (int i=0; i < dim; ++i) {
    var->set_cur(i);
    bool ok = var->get(buffer, 1, bufferSize, 0, 0, 0);
    if (!ok) {
      std::ostringstream msg;
      msg << "Could not read string in variable '" << name << "'.";
      throw std::runtime_error(msg.str());
    } // if
    (*values)[i] = buffer;
    //std::cout << "GROUP: '" << (*values)[i] << "'." << std::endl;
  } // for
  delete[] buffer; buffer = 0;
  delete[] counts; counts = 0;

  PYLITH_METHOD_END;
} // getVar


// End of file 
