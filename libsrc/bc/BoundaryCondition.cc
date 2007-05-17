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

#include "BoundaryCondition.hh" // implementation of object methods

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::BoundaryCondition::BoundaryCondition(void) :
  _id(0),
  _label(""),
  _db(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::BoundaryCondition::~BoundaryCondition(void)
{ // destructor
  _db = 0;
} // destructor

// ----------------------------------------------------------------------
// Set number of degrees of freedom that are constrained at points in field.
void
pylith::bc::BoundaryCondition::setConstraintSizes(
				     const ALE::Obj<real_section_type>& field,
				     const ALE::Obj<ALE::Mesh>& mesh)
{ // setConstraintSizes
} // setConstraintSizes

// ----------------------------------------------------------------------
// Set which degrees of freedom are constrained at points in field.
void
pylith::bc::BoundaryCondition::setConstraints(
				    const ALE::Obj<real_section_type>& field,
				    const ALE::Obj<ALE::Mesh>& mesh)
{ // setConstraints
} // setConstraints

// ----------------------------------------------------------------------
// Set constrained degrees of freedom in field.
void
pylith::bc::BoundaryCondition::integrateJacobian(
				     PetscMat* jacobian,
				     const ALE::Obj<real_section_type>& field,
				     const ALE::Obj<ALE::Mesh>& mesh)
{ // integrateJacobian
} // integrateJacobian

// ----------------------------------------------------------------------
// Set constrained degrees of freedom in field.
void
pylith::bc::BoundaryCondition::integrateResidual(
			      const ALE::Obj<real_section_type>& residual,
			      const ALE::Obj<real_section_type>& fieldT,
			      const ALE::Obj<real_section_type>& fieldTmdt,
			      const ALE::Obj<ALE::Mesh>& mesh)
{ // integrateResidual
} // integrateResidual

// ----------------------------------------------------------------------
// Set values in field.
void
pylith::bc::BoundaryCondition::setField(const double t,
					const ALE::Obj<real_section_type>& field,
					const ALE::Obj<ALE::Mesh>& mesh)
{ // setField
} // setField


// End of file 
