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

#include "FieldBase.hh" // implementation of class methods
#include "pylith/utils/error.hh" // USES std::logic_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::FieldBase::FieldBase(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::FieldBase::~FieldBase(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get string associated with vector field type.
const char*
pylith::topology::FieldBase::vectorFieldString(VectorFieldEnum value)
{ // vectorFieldString
  switch (value) {
  case SCALAR :
    return "scalar";
  case VECTOR :
    return "vector";
  case TENSOR :
    return "tensor";
  case OTHER :
    return "other";
  case MULTI_SCALAR :
    return "multi_scalar";
  case MULTI_VECTOR :
    return "multi_vector";
  case MULTI_TENSOR :
    return "multi_tensor";
  case MULTI_OTHER :
    return "multi_other";
  default :
    assert(0);
    throw std::logic_error("Unknown vector field type in vectorFieldString().");
  } // switch
} // vectorFieldString

// ----------------------------------------------------------------------
// Get string associated with vector field type.
pylith::topology::FieldBase::VectorFieldEnum
pylith::topology::FieldBase::parseVectorFieldString(const char* value)
{ // parseVectorFieldString
  VectorFieldEnum valueEnum = SCALAR;

  if (0 == strcmp(value, "scalar"))
    valueEnum = SCALAR;
  else if (0 == strcmp(value, "vector"))
    valueEnum = VECTOR;
  else if (0 == strcmp(value, "tensor"))
    valueEnum = TENSOR;
  else if (0 == strcmp(value, "other"))
    valueEnum = OTHER;
  else if (0 == strcmp(value, "multi_scalar"))
    valueEnum = MULTI_SCALAR;
  else if (0 == strcmp(value, "multi_vector"))
    valueEnum = MULTI_VECTOR;
  else if (0 == strcmp(value, "multi_tensor"))
    valueEnum = MULTI_TENSOR;
  else if (0 == strcmp(value, "multi_other"))
    valueEnum = MULTI_OTHER;
  else {
    assert(0);
    throw std::logic_error("Unknown vector field string in "
			   "parseVectorFieldString().");
  } // else
  
  return valueEnum;
} // parseVectorFieldString


// End of file 
