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

#include "FaultCohesiveKin.hh" // implementation of object methods

#include "EqKinSrc.hh" // USES EqKinSrc

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveKin::FaultCohesiveKin(void) :
  _eqsrc(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveKin::~FaultCohesiveKin(void)
{ // destructor
  delete _eqsrc; _eqsrc = 0;
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::faults::FaultCohesiveKin::FaultCohesiveKin(const FaultCohesiveKin& f) :
  FaultCohesive(f),
  _eqsrc(0)
{ // copy constructor
  if (0 != f._eqsrc)
    _eqsrc = f._eqsrc->clone();
} // copy constructor

// ----------------------------------------------------------------------
// Set kinematic earthquake source.
void
pylith::faults::FaultCohesiveKin::eqsrc(EqKinSrc* src)
{ // eqsrc
  delete _eqsrc; _eqsrc = (0 != src) ? src->clone() : 0;
} // eqsrc

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveKin::initialize(const ALE::Obj<ALE::Mesh>& mesh,
					     const double_array& upDir)
{ // initialize
  assert(0 != _quadrature);
  
  if (3 != upDir.size())
    throw std::runtime_error("Up direction for fault orientation must be "
			     "a vector with 3 components.");
  
  // Allocate section for orientation information (all vertices)

  // Loop over cells
  //   Compute cell geometry at vertices
  //   Compute weighted orientation of face at vertices (using geometry info)
  //   Update weighted orientations (wt is |J|)

  // Assemble orientation information

  // Loop over vertices
  //   Make orientation information unit magnitude

  // Create list of constraint vertices

  // Only store orientation information at constraint vertices
} // initialize

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term.
void
pylith::faults::FaultCohesiveKin::integrateResidual(
				const ALE::Obj<real_section_type>& residual,
				const ALE::Obj<real_section_type>& disp,
				const ALE::Obj<Mesh>& mesh)
{ // integrateResidual

  // Subtract constraint forces (which are in disp at the constraint
  // DOF) to residual; contributions are at DOF of normal vertices (i and j)

} // integrateResidual

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveKin::integrateJacobian(
				    PetscMat* mat,
				    const ALE::Obj<real_section_type>& dispT,
				    const ALE::Obj<Mesh>& mesh)
{ // integrateJacobian

  // Add constraint information to Jacobian matrix; these are the
  // direction cosines. Entries are associated with vertices ik, jk,
  // ki, and kj.

} // integrateJacobian
  
// ----------------------------------------------------------------------
// Set field.
void
pylith::faults::FaultCohesiveKin::setField(
				     const double t,
				     const ALE::Obj<real_section_type>& disp,
				     const ALE::Obj<Mesh>& mesh)
{ // setField
  typedef std::vector<Mesh::point_type>::const_iterator vert_iterator;

  assert(0 != _eqsrc);

  const ALE::Obj<real_section_type>& slip = _eqsrc->slip(t, _constraintVert);
  assert(!slip.isNull());
  const vert_iterator vBegin = _constraintVert.begin();
  const vert_iterator vEnd = _constraintVert.end();
  for (vert_iterator v_iter=vBegin; v_iter != vEnd; ++v_iter)
    disp->updatePoint(*v_iter, slip->restrictPoint(*v_iter));
} // setField


// End of file 
