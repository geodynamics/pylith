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

#include "pylith/meshio/ExodusII.hh" // implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES journal macros
#include "pylith/utils/Exceptions.hh" // USES Exception

#include "petsc.h" // USES MPI_Comm

// :KLUDGE: Prevent NetCDF from definining MPI types
#define MPI_INCLUDED
#include <netcdf.h> // USES netcdf

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::ExodusII::ExodusII(void) :
    _filename(""),
    _file(0) { // constructor
} // constructor


// ----------------------------------------------------------------------
// Constructor with filename
pylith::meshio::ExodusII::ExodusII(const char* filename) :
    _filename(filename),
    _file(0) { // constructor
    open();
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::ExodusII::~ExodusII(void) { // destructor
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate data structures.
void
pylith::meshio::ExodusII::deallocate(void) { // deallocate
    PYLITH_METHOD_BEGIN;

    close();

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Set filename.
void
pylith::meshio::ExodusII::filename(const char* name) { // filename
    _filename = name;
} // filename


// ----------------------------------------------------------------------
// Get filename.
const char*
pylith::meshio::ExodusII::filename(void) const { // filename
    return _filename.c_str();
} // filename


// ----------------------------------------------------------------------
// Open file.
void
pylith::meshio::ExodusII::open(void) { // open
    PYLITH_METHOD_BEGIN;

    close();

    const int err = nc_open(_filename.c_str(), NC_NOWRITE, &_file);
    if (err != NC_NOERR) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Could not open ExodusII file '" << _filename << "' for reading.");
    } // if

    PYLITH_METHOD_END;
} // open


// ----------------------------------------------------------------------
// Close file.
void
pylith::meshio::ExodusII::close(void) { // close
    PYLITH_METHOD_BEGIN;

    if (_file) {
        int err = nc_close(_file);_file = 0;
        if (err != NC_NOERR) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                         "Could not close ExodusII file '" << _filename << "'.");
        } // if
    } // if

    PYLITH_METHOD_END;
} // close


// ----------------------------------------------------------------------
// Check if Cubit Exodus file constains dimension.
bool
pylith::meshio::ExodusII::hasDim(const char* name,
                                 int* id) const { // hasDim
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
pylith::meshio::ExodusII::hasAttr(const char* name,
                                  const char* varName,
                                  int* id) const {
    PYLITH_METHOD_BEGIN;

    assert(_file);

    bool found = false;
    if (id) {
        *id = -1; // Not found
    } // if

    int varId = -1;
    if (varName && !hasVar(varName, &varId)) {
        PYLITH_METHOD_RETURN(found);
    } // if

    int attrId = -1;
    if (NC_NOERR == nc_inq_attid(_file, varId, name, &attrId)) {
        found = true;
        if (id) {
            *id = attrId;
        } // if
    } // if

    PYLITH_METHOD_RETURN(found);
} // hasAttr


// ----------------------------------------------------------------------
// Check if Cubit Exodus file constains variable.
bool
pylith::meshio::ExodusII::hasVar(const char* name,
                                 int* id) const {
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
    } else if (id) {
        *id = -1;
    } // if/else

    PYLITH_METHOD_RETURN(found);
} // hasVar


// ----------------------------------------------------------------------
// Get value for dimension.
int
pylith::meshio::ExodusII::getDim(const char* name) const { // getDim
    PYLITH_METHOD_BEGIN;

    assert(_file);

    int id = -1;
    if (!hasDim(name, &id)) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Missing dimension '" << name << "'.");
    } // if

    size_t dimSize = 0;
    int err = nc_inq_dimlen(_file, id, &dimSize);
    if (err != NC_NOERR) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Could not get dimension '" << name << "'.");
    } // if

    PYLITH_METHOD_RETURN(dimSize);
} // getDim


// ----------------------------------------------------------------------
// Get string attribute.
void
pylith::meshio::ExodusII::getAttr(std::string* value,
                                  const char* varName,
                                  const char* name) const {
    PYLITH_METHOD_BEGIN;

    assert(_file);
    assert(value);

    int varId = NC_GLOBAL;
    if (varName && !hasVar(varName, &varId)) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Could not find variable '" << varName << "' to get string attribute '" << name << "'.");
    } // if
    assert(!varName || varId >= 0);

    int attrId = -1;
    if (!hasAttr(name, varName, &attrId)) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Missing string attribute '" << name << "'.");
    } // if

    int err = NC_NOERR;
    size_t stringLen = 0;
    err = nc_inq_attlen(_file, varId, name, &stringLen);
    if (NC_NOERR != err) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Could not get length of string attribute '" << name << "' for variable '" << varName << "'.");
    } // if
    char* cValue = stringLen > 0 ? new char[stringLen+1] : nullptr;
    if (cValue) {
        err = nc_get_att_text(_file, varId, name, cValue);
        if (err != NC_NOERR) {
            PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                         "Could not read string attribute '" << name << "' for variable '" << varName << "'.");
        } // if
        *value = cValue;
        delete[] cValue;cValue = nullptr;
    } // if

    PYLITH_METHOD_END;
} // getAttr


// ----------------------------------------------------------------------
// Get values for variable as an array of PylithScalars.
void
pylith::meshio::ExodusII::getVar(PylithScalar* values,
                                 int* dims,
                                 int ndims,
                                 const char* name) const { // getVar
    PYLITH_METHOD_BEGIN;

    assert(_file);
    assert(values);

    int vid = -1;
    if (!hasVar(name, &vid)) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Missing real variable '" << name << "'.");
    } // if

    int vndims = 0;
    int err = nc_inq_varndims(_file, vid, &vndims);
    if (ndims != vndims) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Expecting " << ndims << " dimensions for variable '" << name
                                  << "' but variable only has " << vndims << " dimensions.");
    } // if

    int* dimIds = (ndims > 0) ? new int[ndims] : 0;
    err = nc_inq_vardimid(_file, vid, dimIds);
    if (err != NC_NOERR) {
        delete[] dimIds;dimIds = 0;
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Could not get dimensions for variable '" << name << "'.");
    } // if

    for (int iDim = 0; iDim < ndims; ++iDim) {
        size_t dimSize = 0;
        err = nc_inq_dimlen(_file, dimIds[iDim], &dimSize);
        if (err != NC_NOERR) {
            delete[] dimIds;dimIds = 0;
            PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                         "Could not get dimension '" << iDim << "' for variable '" << name << "'.");
        } // if
        if (size_t(dims[iDim]) != dimSize) {
            delete[] dimIds;dimIds = 0;
            PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                         "Expecting dimension " << iDim << " of variable '" << name
                                                << "' to be " << dims[iDim] << ", but dimension is " << dimSize
                                                << ".");
        } // if
    } // for
    delete[] dimIds;dimIds = 0;

    if (sizeof(PylithScalar) == sizeof(double)) {
        err = nc_get_var_double(_file, vid, values);
    } else {
        PYLITH_ERROR(pylith::InternalLogicError, pylith::journal::logic,
                     "Unknown size of PylithScalar in ExodusII::getVar().");
    } // if/else
    if (err != NC_NOERR) {
        PYLITH_ERROR(pylith::InternalLogicError, pylith::journal::logic,
                     "Could not get values for variable '" << name << ".");
    } // if

    PYLITH_METHOD_END;
} // getVar


// ----------------------------------------------------------------------
// Get values for variable as an array of ints.
void
pylith::meshio::ExodusII::getVar(int* values,
                                 int* dims,
                                 int ndims,
                                 const char* name) const { // getVar
    PYLITH_METHOD_BEGIN;

    assert(_file);
    assert(values);

    int vid = -1;
    if (!hasVar(name, &vid)) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Missing integer variable '" << name << "'.");
    } // if

    int vndims = 0;
    int err = nc_inq_varndims(_file, vid, &vndims);
    if (ndims != vndims) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Expecting " << ndims << " dimensions for variable '" << name
                                  << "' but variable only has " << vndims << " dimensions.");
    } // if

    int* dimIds = (ndims > 0) ? new int[ndims] : 0;
    err = nc_inq_vardimid(_file, vid, dimIds);
    if (err != NC_NOERR) {
        delete[] dimIds;dimIds = 0;
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Could not get dimensions for variable '" << name << "'.");
    } // if

    for (int iDim = 0; iDim < ndims; ++iDim) {
        size_t dimSize = 0;
        err = nc_inq_dimlen(_file, dimIds[iDim], &dimSize);
        if (err != NC_NOERR) {
            delete[] dimIds;dimIds = 0;
            PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                         "Could not get dimension '" << iDim << "' for variable '" << name << "'.");
        } // if
        if (size_t(dims[iDim]) != dimSize) {
            delete[] dimIds;dimIds = 0;
            PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                         "Expecting dimension " << iDim << " of variable '" << name
                                                << "' to be " << dims[iDim] << ", but dimension is " << dimSize
                                                << ".");
        } // if
    } // for
    delete[] dimIds;dimIds = 0;

    err = nc_get_var_int(_file, vid, values);
    if (err != NC_NOERR) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Could not get values for variable '" << name << ".");
    } // if

    PYLITH_METHOD_END;
} // getVar


// ----------------------------------------------------------------------
// Get values for variable as an array of strings.
void
pylith::meshio::ExodusII::getVar(string_vector* values,
                                 int dim,
                                 const char* name) const { // getVar
    PYLITH_METHOD_BEGIN;

    assert(_file);
    assert(values);

    int vid = -1;
    if (!hasVar(name, &vid)) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Missing string variable '" << name << "'.");
    } // if

    int vndims = 0;
    const int ndims = 2;
    int err = nc_inq_varndims(_file, vid, &vndims);
    if (ndims != vndims) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Expecting " << ndims << " dimensions for variable '" << name
                                  << "' but variable only has " << vndims << " dimensions.");
    } // if

    int dimIds[ndims];
    err = nc_inq_vardimid(_file, vid, dimIds);
    if (err != NC_NOERR) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Could not get dimensions for variable '" << name << "'.");
    } // if

    int iDim = 0;
    size_t dimSize = 0;
    err = nc_inq_dimlen(_file, dimIds[0], &dimSize);
    if (err != NC_NOERR) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Could not get dimension '" << iDim << "' for variable '" << name << "'.");
    } // if
    if (size_t(dim) != dimSize) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Expecting dimension " << iDim << " of variable '" << name
                                            << "' to be " << dim << ", but dimension is " << dimSize
                                            << ".");
    } // if

    iDim = 1;
    err = nc_inq_dimlen(_file, dimIds[1], &dimSize);
    if (err != NC_NOERR) {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                     "Could not get dimension '" << iDim << "' for variable '" << name << "'.");
    } // if
    const size_t bufferSize = dimSize;
    char* buffer = (bufferSize > 0) ? new char[bufferSize] : 0;
    size_t indices[2] = { 0, 0 };
    size_t chunk[2] = { 1, bufferSize };
    for (int i = 0; i < dim; ++i) {
        indices[0] = i;
        err = nc_get_vara_text(_file, vid, indices, chunk, buffer);
        if (err != NC_NOERR) {
            delete[] buffer;buffer = 0;
            PYLITH_ERROR(pylith::IOError, pylith::journal::user_input,
                         "Could not read string in variable '" << name << "'.");
        } // if
        (*values)[i] = buffer;
        // std::cout << "GROUP: '" << (*values)[i] << "'." << std::endl;
    } // for
    delete[] buffer;buffer = nullptr;

    PYLITH_METHOD_END;
} // getVar


// End of file
