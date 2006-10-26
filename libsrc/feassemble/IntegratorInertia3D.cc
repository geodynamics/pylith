// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "IntegratorInertia3D.hh" // implementation of class methods

#include "petscmat.h" // USES PetscMat

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorInertia3D::IntegratorInertia3D(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorInertia3D::~IntegratorInertia3D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::IntegratorInertia3D::IntegratorInertia3D(const IntegratorInertia3D& i) :
  Integrator(i)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Integrate inertial term for 3-D finite elements.
void
pylith::feassemble::IntegratorInertia3D::integrateAction(const ALE::Obj<ALE::Mesh::real_section_type>& fieldOut,
		 const ALE::Obj<ALE::Mesh::real_section_type>& fieldIn,
		 const ALE::Obj<ALE::Mesh::real_section_type>& coordinates)
{ // integrateAction
} // integrateAction

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::IntegratorInertia3D::integrate(PetscMat* mat,
		    const ALE::Obj<ALE::Mesh::real_section_type>& coordinates)
{ // integrate
} // integrate


// End of file 
