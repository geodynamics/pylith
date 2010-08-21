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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "SlipTimeFn.hh" // implementation of object methods

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Field.hh" // USES Field

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::SlipTimeFn::SlipTimeFn(void) :
  _parameters(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::SlipTimeFn::~SlipTimeFn(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::SlipTimeFn::deallocate(void)
{ // deallocate
  delete _parameters; _parameters = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Get parameter fields.
const pylith::topology::Fields<pylith::topology::Field<pylith::topology::SubMesh> >*
pylith::faults::SlipTimeFn::parameterFields(void) const
{ // parameterFields
  return _parameters;
} // parameterFields


// End of file 
