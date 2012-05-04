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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "Nucleator.hh" // implementation of object methods

#include "TractPerturbation.hh" // USES SlipTimeFn

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::Nucleator::Nucleator(void) :
  _originTime(0.0),
  _tractfn(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::Nucleator::~Nucleator(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::Nucleator::deallocate(void)
{ // deallocate
  _tractfn = 0; // :TODO: Use shared pointer.
} // deallocate
  
// ----------------------------------------------------------------------
// Set origin time for earthquake source.
void
pylith::faults::Nucleator::originTime(const PylithScalar value)
{ // originTime
  _originTime = value;
} // originTime

// ----------------------------------------------------------------------
// Get origin time for earthquake source.
PylithScalar
pylith::faults::Nucleator::originTime(void) const
{ // originTime
  return _originTime;
} // originTime

// ----------------------------------------------------------------------
// Set slip time function.
void
pylith::faults::Nucleator::perturbationFn(TractPerturbation* tractfn)
{ // perturbationFn
  _tractfn = tractfn; // :TODO: Use shared pointer.
} // perturbationFn

// ----------------------------------------------------------------------
// Initialize slip time function.
void
pylith::faults::Nucleator::initialize(
			   const topology::SubMesh& faultMesh,
			   const spatialdata::units::Nondimensional& normalizer)
{ // initialize
  // :TODO: Normalize origin time in Python?
  normalizer.nondimensionalize(&_originTime, 1, normalizer.timeScale());
  assert(_tractfn);
  _tractfn->initialize(faultMesh, normalizer);
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
void
pylith::faults::Nucleator::traction(
			   topology::Field<topology::SubMesh>* const tractionField,
			   const PylithScalar t)
{ // slip
  assert(_tractfn);
  _tractfn->traction(tractionField, t-_originTime);
} // slip

// ----------------------------------------------------------------------
// Get amplitude of spatial variation of traction.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::Nucleator::tractionAmp(void) const
{ // finalSlip
  assert(_tractfn);
  return _tractfn->amplitude();
} // tractionAmp


// End of file 
