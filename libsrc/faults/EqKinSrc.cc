// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "EqKinSrc.hh" // implementation of object methods

#include "SlipTimeFn.hh" // USES SlipTimeFn

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <assert.h> // USES assert()

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
  _slipfn = 0; // Don't manage memory for slip fn
} // destructor

// ----------------------------------------------------------------------
// Set origin time for earthquake source.
void
pylith::faults::EqKinSrc::originTime(const double value)
{ // originTime
  _originTime = value;
} // originTime

// ----------------------------------------------------------------------
// Get origin time for earthquake source.
double
pylith::faults::EqKinSrc::originTime(void) const
{ // originTime
  return _originTime;
} // originTime

// ----------------------------------------------------------------------
// Set slip time function.
void
pylith::faults::EqKinSrc::slipfn(SlipTimeFn* slipfn)
{ // slipfn
  _slipfn = slipfn; // Don't manage memory for slip fn
} // slipfn

// ----------------------------------------------------------------------
// Initialize slip time function.
void
pylith::faults::EqKinSrc::initialize(
			 const ALE::Obj<Mesh>& faultMesh,
			 const spatialdata::geocoords::CoordSys* cs,
			 const spatialdata::units::Nondimensional& normalizer)
{ // initialize
  normalizer.nondimensionalize(&_originTime, 1, normalizer.timeScale());
  assert(0 != _slipfn);
  _slipfn->initialize(faultMesh, cs, normalizer, _originTime);
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
void
pylith::faults::EqKinSrc::slip(const ALE::Obj<pylith::real_section_type>& slipField,
			       const double t,
			       const ALE::Obj<Mesh>& faultMesh)
{ // slip
  assert(0 != _slipfn);
  _slipfn->slip(slipField, t, faultMesh);
} // slip

// ----------------------------------------------------------------------
// Get slip increment on fault surface from time t0 to 1.
void
pylith::faults::EqKinSrc::slipIncr(const ALE::Obj<pylith::real_section_type>& slipField,
				   const double t0,
				   const double t1,
				   const ALE::Obj<Mesh>& faultMesh)
{ // slip
  assert(0 != _slipfn);
  _slipfn->slipIncr(slipField, t0, t1, faultMesh);
} // slip

// ----------------------------------------------------------------------
// Get final slip.
ALE::Obj<pylith::real_section_type>
pylith::faults::EqKinSrc::finalSlip(void)
{ // finalSlip
  assert(0 != _slipfn);
  return _slipfn->finalSlip();
} // finalSlip

// ----------------------------------------------------------------------
// Get time when slip begins at each point.
ALE::Obj<pylith::real_section_type>
pylith::faults::EqKinSrc::slipTime(void)
{ // slipTime
  assert(0 != _slipfn);
  return _slipfn->slipTime();
} // slipTime


// End of file 
