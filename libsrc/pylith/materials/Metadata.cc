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

#include "Metadata.hh" // implementation of class methods

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error
#include <iostream> // USES std::cerr

// ----------------------------------------------------------------------
// Constructor.
pylith::materials::Metadata::Metadata(const ParamDescription* props,
				      const int numProps,
				      const char* dbProps[],
				      const int numDBProps,
				      const ParamDescription* vars,
				      const int numVars,
				      const char* dbVars[],
				      const int numDBVars) :
  _properties(0),
  _stateVars(0),
  _dbProperties(dbProps),
  _dbStateVars(dbVars),
  _numProperties(numProps),
  _numStateVars(numVars),
  _numDBProperties(numDBProps),
  _numDBStateVars(numDBVars)
{ // constructor
  if (numProps > 0) {
    assert(props);
    _properties = new ParamDescription[numProps];
    for (int i=0; i < numProps; ++i)
      _properties[i] = props[i];
  } // if

  if (numVars > 0) {
    assert(vars);
    _stateVars = new ParamDescription[numVars];
    for (int i=0; i < numVars; ++i)
      _stateVars[i] = vars[i];
  } // if
} // constructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::materials::Metadata::Metadata(const Metadata& m) :
  _properties(0),
  _stateVars(0),
  _dbProperties(m._dbProperties),
  _dbStateVars(m._dbStateVars),
  _numProperties(m._numProperties),
  _numStateVars(m._numStateVars),
  _numDBProperties(m._numDBProperties),
  _numDBStateVars(m._numDBStateVars)
{ // copy constructor
  if (m._properties) {
    assert(_numProperties > 0);
    _properties = new ParamDescription[_numProperties];
    for (int i=0; i < _numProperties; ++i)
      _properties[i] = m._properties[i];
  } // if

  if (m._stateVars) {
    assert(_numStateVars > 0);
    _stateVars = new ParamDescription[_numStateVars];
    for (int i=0; i < _numStateVars; ++i)
      _stateVars[i] = m._stateVars[i];
  } // if
} // copy constructor

// ----------------------------------------------------------------------
// Default destructor
pylith::materials::Metadata::~Metadata(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::Metadata::deallocate(void)
{ // deallocate
  delete[] _properties; _properties = 0;
  delete[] _stateVars; _stateVars = 0;
} // deallocate
  

// End of file
