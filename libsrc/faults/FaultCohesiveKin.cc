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
  
  
  
} // initialize

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term.
void
pylith::faults::FaultCohesiveKin::integrateResidual(
				const ALE::Obj<real_section_type>& residual,
				const ALE::Obj<real_section_type>& disp,
				const ALE::Obj<Mesh>& mesh)
{ // integrateResidual
} // integrateResidual

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveKin::integrateJacobian(
				    PetscMat* mat,
				    const ALE::Obj<real_section_type>& dispT,
				    const ALE::Obj<Mesh>& mesh)
{ // integrateJacobian
} // integrateJacobian
  
// ----------------------------------------------------------------------
// Set field.
void
pylith::faults::FaultCohesiveKin::setField(
				     const ALE::Obj<real_section_type>& disp,
				     const ALE::Obj<Mesh>& mesh)
{ // setField
} // setField


// End of file 
