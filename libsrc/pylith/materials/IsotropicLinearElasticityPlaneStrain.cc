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

#include "IsotropicLinearElasticityPlaneStrain.hh" // implementation of object methods

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearElasticityPlaneStrain::IsotropicLinearElasticityPlaneStrain(void) :
  Material(2)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearElasticityPlaneStrain::~IsotropicLinearElasticityPlaneStrain(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Include inertia?
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::useInertia(const bool value)
{ // useInertia
  _useInertia = value;
} // useInertia

// ----------------------------------------------------------------------
// Include body force?
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::useBodyForce(const bool value)
{ // useBodyForce
  _useBodyForce = value;
} // useBodyForce

// ----------------------------------------------------------------------
// Preinitialize material. Set names/sizes of auxiliary fields.
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::preinitialize(void)
{ // preinitialize
  // Set db values
  // :TODO: ADD STUFF HERE, values depend on _useInertia, _useBodyForce
  
  // Set auxiliary fields
  // :TODO: ADD STUFF HERE, values depend on _useInertia, _useBodyForce
} // preinitialize

// ----------------------------------------------------------------------
// Set residual and Jacobian kernels.
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernels(const PetscDS prob) const
{ // _setFEKernels
  const PetscInt disp = 0; // :KLUDGE: These need to come from the problem!!!!!
  const PetscInt vel = 1;

  // :TODO: Select kernels based on _useIneria and _useBodyForce

  // Residual kernels
  const PetscPointResFunc f0_disp = NULL;
  const PetscPointResFunc f1_disp = NULL;
  err = PetscDSSetResidual(prob, disp, f0_disp, f1_disp);PYLITH_CHECK_ERROR(err);

  const PetscPointResFunc f0_vel = NULL;
  const PetscPointResFunc f1_vel = NULL;
  err = PetscDSSetResidual(prob, vel, f0_vel, f1_vel);PYLITH_CHECK_ERROR(err);

  // Jacobian kernels
  const PetscPointJacFunc g0_dispdisp = NULL;
  const PetscPointJacFunc g1_dispdisp = NULL;
  const PetscPointJacFunc g2_dispdisp = NULL;
  const PetscPointJacFunc g3_dispdisp = NULL;
  const PetscPointJacFunc g0_dispvel = NULL;
  const PetscPointJacFunc g1_dispvel = NULL;
  const PetscPointJacFunc g2_dispvel = NULL;
  const PetscPointJacFunc g3_dispvel = NULL;
  const PetscPointJacFunc g0_veldisp = NULL;
  const PetscPointJacFunc g1_veldisp = NULL;
  const PetscPointJacFunc g2_veldisp = NULL;
  const PetscPointJacFunc g3_veldisp = NULL;
  const PetscPointJacFunc g0_velvel = NULL;
  const PetscPointJacFunc g1_velvel = NULL;
  const PetscPointJacFunc g2_velvel = NULL;
  const PetscPointJacFunc g3_velvel = NULL;
  err = PetscDSSetJacobian(prob, disp, disp, g0_dispdisp, g1_dispdisp, g2_dispdisp, g3_dispdisp);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, disp, vel, g0_dispvel, g1_dispvel, g2_dispvel, g3_dispvel);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, vel, disp, g0_veldisp, g1_veldisp, g2_veldisp, g3_veldisp);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, vel, vel, g0_velvel, g1_velvel, g2_velvel, g3_velvel);PYLITH_CHECK_ERROR(err);

} // _setFEKernels

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_dbToAuxFields(PylithScalar const auxValues[],
									const int numAuxValues,
									const scalar_array& dbValues) const
{ // _dbToAuxFields
  // :TODO: ADD STUFF HERE, values depend on _useInertia, _useBodyForce
} // _dbToaAuxFields
  
// ----------------------------------------------------------------------
// Nondimensionalize auxiliary fields. Nondimensionalization is done
// in place (no copy).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_nondimAuxFields(PylithScalar const values[],
									  const int nvalues) const
{ // _nondimAuxFields
  // :TODO: ADD STUFF HERE, use scales in _auxFields
  // Can we move this to Material?
} // _nondimAuxFields
  
// ----------------------------------------------------------------------
// Dimensionalize auxiliary fields. Dimensionalization is done in
// place (no copy).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_dimAuxFields(PylithScalar const values[],
								       const int nvalues) const
{ // _dimAuxFields
  // :TODO: ADD STUFF HERE, use scales in _auxFields
  // Can we move this to Material?
} // _dimAuxFields


// End of file 
