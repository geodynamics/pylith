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

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::EqKinSrc::EqKinSrc(void) :
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
			      const ALE::Obj<Mesh>& mesh,
			      const ALE::Obj<Mesh>& faultMesh,
			      const std::set<Mesh::point_type>& vertices,
			      const spatialdata::geocoords::CoordSys* cs)
{ // initialize
  assert(0 != _slipfn);
  _slipfn->initialize(mesh, faultMesh, vertices, cs);
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
const ALE::Obj<pylith::real_section_type>&
pylith::faults::EqKinSrc::slip(const double t,
			       const std::set<Mesh::point_type>& vertices)
{ // slip
  assert(0 != _slipfn);
  return _slipfn->slip(t, vertices);
} // slip

// ----------------------------------------------------------------------
// Get slip increment on fault surface from time t0 to 1.
const ALE::Obj<pylith::real_section_type>&
pylith::faults::EqKinSrc::slipIncr(const double t0,
				   const double t1,
				   const std::set<Mesh::point_type>& vertices)
{ // slip
  assert(0 != _slipfn);
  return _slipfn->slipIncr(t0, t1, vertices);
} // slip


// End of file 
