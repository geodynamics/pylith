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

#include "EqKinSrc.hh" // implementation of object methods

#include "SlipTimeFn.hh" // USES SlipTimeFn

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::EqKinSrc::EqKinSrc(void) :
  _originTime(0.0),
  _slipfn(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::EqKinSrc::~EqKinSrc(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::EqKinSrc::deallocate(void)
{ // deallocate
  _slipfn = 0; // :TODO: Use shared pointer.
} // deallocate
  
// ----------------------------------------------------------------------
// Set origin time for earthquake source.
void
pylith::faults::EqKinSrc::originTime(const PylithScalar value)
{ // originTime
  _originTime = value;
} // originTime

// ----------------------------------------------------------------------
// Get origin time for earthquake source.
PylithScalar
pylith::faults::EqKinSrc::originTime(void) const
{ // originTime
  return _originTime;
} // originTime

// ----------------------------------------------------------------------
// Set slip time function.
void
pylith::faults::EqKinSrc::slipfn(SlipTimeFn* slipfn)
{ // slipfn
  _slipfn = slipfn; // :TODO: Use shared pointer.
} // slipfn

// ----------------------------------------------------------------------
// Initialize slip time function.
void
pylith::faults::EqKinSrc::initialize(const topology::Mesh& faultMesh,
				     const spatialdata::units::Nondimensional& normalizer)
{ // initialize
  PYLITH_METHOD_BEGIN;

  // :TODO: Normalize slip time in Python?
  normalizer.nondimensionalize(&_originTime, 1, normalizer.timeScale());
  assert(_slipfn);
  _slipfn->initialize(faultMesh, normalizer, _originTime);

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
void
pylith::faults::EqKinSrc::slip(topology::Field* const slipField,
			       const PylithScalar t)
{ // slip
  PYLITH_METHOD_BEGIN;

  assert(_slipfn);
  _slipfn->slip(slipField, t);

  PYLITH_METHOD_END;
} // slip

// ----------------------------------------------------------------------
// Get final slip.
const pylith::topology::Field&
pylith::faults::EqKinSrc::finalSlip(void) const
{ // finalSlip
  PYLITH_METHOD_BEGIN;

  assert(_slipfn);
  PYLITH_METHOD_RETURN(_slipfn->finalSlip());
} // finalSlip

// ----------------------------------------------------------------------
// Get time when slip begins at each point.
const pylith::topology::Field&
pylith::faults::EqKinSrc::slipTime(void) const
{ // slipTime
  PYLITH_METHOD_BEGIN;

  assert(_slipfn);
  PYLITH_METHOD_RETURN(_slipfn->slipTime());
} // slipTime


// End of file 
