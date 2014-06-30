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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "SolutionFields.hh" // implementation of class methods

#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::SolutionFields::SolutionFields(const Mesh& mesh) :
  Fields(mesh),
  _solutionName("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::SolutionFields::~SolutionFields(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::SolutionFields::deallocate(void)
{ // deallocate
  Fields::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Set name of solution field.
void
pylith::topology::SolutionFields::solutionName(const char* name)
{ // solutionName
  map_type::const_iterator iter = _fields.find(name);
  if (iter == _fields.end()) {
    std::ostringstream msg;
    msg << "Cannot use unknown field '" << name << "' when setting name of solution field.";
    throw std::runtime_error(msg.str());
  } // if
  _solutionName = name;
} // solutionName

// ----------------------------------------------------------------------
// Get solution field.
const pylith::topology::Field&
pylith::topology::SolutionFields::solution(void) const
{ // solution
  if (_solutionName == "")
    throw std::runtime_error("Cannot retrieve solution. Name of solution " \
			     "field has not been specified.");
  return get(_solutionName.c_str());
} // solution

// ----------------------------------------------------------------------
// Get solution field.
pylith::topology::Field&
pylith::topology::SolutionFields::solution(void)
{ // solution
  if (_solutionName == "")
    throw std::runtime_error("Cannot retrieve solution. Name of solution " \
			     "field has not been specified.");
  return get(_solutionName.c_str());
} // solution


// End of file 
