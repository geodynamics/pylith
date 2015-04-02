// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "Constraint.hh" // implementation of object methods

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::Constraint::Constraint(void) :
  _normalizer(new spatialdata::units::Nondimensional)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::feassemble::Constraint::~Constraint(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::Constraint::deallocate(void)
{ // deallocate
  delete _normalizer; _normalizer = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::feassemble::Constraint::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
  if (!_normalizer)
    _normalizer = new spatialdata::units::Nondimensional(dim);
  else
    *_normalizer = dim;
} // normalizer


// End of file 
