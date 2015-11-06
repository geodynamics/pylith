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

#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/materials/Query.hh" // USES Query
extern "C" {
#include "pylith/fekernels/elasticity.h" // USES fekernels
}

#include "petscds.h"

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearElasticityPlaneStrain::IsotropicLinearElasticityPlaneStrain(void) :
  MaterialNew(2),
  _useInertia(false),
  _useBodyForce(false)
{ // constructor
  _isJacobianSymmetric = true;
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
  if (_useInertia) {
    _isJacobianSymmetric = false; 
  } else {
    _isJacobianSymmetric = true;
  } // if/else
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
pylith::materials::IsotropicLinearElasticityPlaneStrain::_auxFieldsSetup(void)
{ // preinitialize
  PYLITH_METHOD_BEGIN;

  // Set subfields in auxiliary fields.
  assert(_normalizer);
  const PylithReal densityScale = _normalizer->densityScale();
  const PylithReal pressureScale = _normalizer->pressureScale();
  const PylithReal lengthScale = _normalizer->lengthScale();
  const PylithReal timeScale = _normalizer->timeScale();
  const PylithReal forceScale = densityScale * lengthScale / (timeScale * timeScale);

  // :ATTENTION: The order here must match the order of the auxiliary fields in the FE kernels.
  _auxFields->subfieldAdd("density", 1, topology::Field::SCALAR, this->discretization("density"), densityScale);
  _auxFields->subfieldAdd("mu", 1, topology::Field::SCALAR, this->discretization("mu"), pressureScale);
  _auxFields->subfieldAdd("lambda", 1, topology::Field::SCALAR, this->discretization("lambda"), pressureScale);
  if (_useBodyForce) {
    _auxFields->subfieldAdd("body force", dimension(), topology::Field::VECTOR, this->discretization("body force"), forceScale);
  } // if

  // Order does not matter.
  _auxFieldsQuery->queryFn("density", pylith::materials::Query::dbQueryDensity2D);
  _auxFieldsQuery->queryFn("mu", pylith::materials::Query::dbQueryMu2D);
  _auxFieldsQuery->queryFn("lambda", pylith::materials::Query::dbQueryLambda2D);
  if (_useBodyForce) {
    _auxFieldsQuery->queryFn("body force", pylith::materials::Query::dbQueryBodyForce2D);
  } // if

  PYLITH_METHOD_END;
} // preinitialize

// ----------------------------------------------------------------------
// Set kernels for RHS residual G(t,u).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernelsRHSResidual(const topology::Field& solution) const
{ // _setFEKernelsRHSResidual
  PYLITH_METHOD_BEGIN;

  const PetscDM dm = solution.dmMesh();assert(dm);
  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(dm, &prob);PYLITH_CHECK_ERROR(err);

  // Displacement
  const PetscInt i_disp = solution.subfieldInfo("displacement").index;
  const PetscPointFunc g0_disp = (_useBodyForce) ? pylith_fekernels_g0_ElasticityBodyForce : NULL;
  const PetscPointFunc g1_disp = pylith_fekernels_g1_IsotropicLinearElasticityPlaneStrain;
  err = PetscDSSetResidual(prob, i_disp, g0_disp, g1_disp);PYLITH_CHECK_ERROR(err);

  // Velocity
  const PetscInt i_vel = solution.subfieldInfo("velocity").index;
  const PetscPointFunc g0_vel = pylith_fekernels_g0_DispVel;
  const PetscPointFunc g1_vel = NULL;
  err = PetscDSSetResidual(prob, i_vel, g0_vel, g1_vel);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _setFEKernelsRHSResidual


// ----------------------------------------------------------------------
// Set kernels for RHS Jacobian G(t,u).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernelsRHSJacobian(const topology::Field& solution) const
{ // _setFEKernelsRHSJacobian
  PYLITH_METHOD_BEGIN;

  const PetscInt i_disp = solution.subfieldInfo("displacement").index;
  const PetscInt i_vel = solution.subfieldInfo("velocity").index;

  const PetscDM dm = solution.dmMesh();assert(dm);
  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(dm, &prob);PYLITH_CHECK_ERROR(err);

  // Jacobian kernels
  const PetscPointJac hg0_dispdisp = NULL;
  const PetscPointJac hg1_dispdisp = NULL;
  const PetscPointJac hg2_dispdisp = NULL;
  const PetscPointJac hg3_dispdisp = pylith_fekernels_hg3_dispdisp_IsotropicLinearElasticityPlaneStrain;
  
  const PetscPointJac hg0_dispvel = NULL;
  const PetscPointJac hg1_dispvel = NULL;
  const PetscPointJac hg2_dispvel = NULL;
  const PetscPointJac hg3_dispvel = NULL;

  const PetscPointJac hg0_veldisp = NULL;
  const PetscPointJac hg1_veldisp = NULL;
  const PetscPointJac hg2_veldisp = NULL;
  const PetscPointJac hg3_veldisp = NULL;

  const PetscPointJac hg0_velvel = pylith_fekernels_hg0_velvel_DispVel;
  const PetscPointJac hg1_velvel = NULL;
  const PetscPointJac hg2_velvel = NULL;
  const PetscPointJac hg3_velvel = NULL;
    
  err = PetscDSSetJacobian(prob, i_disp, i_disp, hg0_dispdisp, hg1_dispdisp, hg2_dispdisp, hg3_dispdisp);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_disp, i_vel,  hg0_dispvel,  hg1_dispvel,  hg2_dispvel,  hg3_dispvel);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_disp, hg0_veldisp,  hg1_veldisp,  hg2_veldisp,  hg3_veldisp);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_vel,  hg0_velvel,   hg1_velvel,   hg2_velvel,   hg3_velvel);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobian


// ----------------------------------------------------------------------
// Set kernels for LHS residual F(t,u,\dot{u}).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernelsLHSResidual(const topology::Field& solution) const
{ // _setFEKernelsLHSResidual
  PYLITH_METHOD_BEGIN;

  const PetscDM dm = solution.dmMesh();assert(dm);
  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(dm, &prob);PYLITH_CHECK_ERROR(err);

  // Displacement
  const PetscInt i_disp = solution.subfieldInfo("displacement").index;
  const PetscPointFunc f0_disp = (_useInteria) ? pylith_fekernels_f0_ElasticityInertia : NULL;
  const PetscPointFunc f1_disp = NULL;
  err = PetscDSSetResidual(prob, i_disp, f0_disp, f1_disp);PYLITH_CHECK_ERROR(err);

  // Velocity
  const PetscInt i_vel = solution.subfieldInfo("velocity").index;
  const PetscPointFunc f0_vel = pylith_fekernels_f0_DispVel;
  const PetscPointFunc f1_vel = NULL;
  err = PetscDSSetResidual(prob, i_vel, f0_vel, f1_vel);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _setFEKernelsLHSResidual


// ----------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,u,\dot{u}).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernelsLHSJacobianImplicit(const topology::Field& solution) const
{ // _setFEKernelsLHSJacobianImplicit
  PYLITH_METHOD_BEGIN;

  const PetscInt i_disp = solution.subfieldInfo("displacement").index;
  const PetscInt i_vel = solution.subfieldInfo("velocity").index;

  const PetscDM dm = solution.dmMesh();assert(dm);
  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(dm, &prob);PYLITH_CHECK_ERROR(err);

  // Jacobian kernels
  const PetscPointJac hf0_dispdisp = NULL;
  const PetscPointJac hf1_dispdisp = NULL;
  const PetscPointJac hf2_dispdisp = NULL;
  const PetscPointJac hf3_dispdisp = NULL;
  
  const PetscPointJac hf0_dispvel = (_useInertia) ? pylith_fekernels_hf0_dispdisp_ElasticityInertiaIm : NULL;
  const PetscPointJac hf1_dispvel = NULL;
  const PetscPointJac hf2_dispvel = NULL;
  const PetscPointJac hf3_dispvel = NULL;

  const PetscPointJac hf0_veldisp = pylith_fekernels_hf0_veldisp_DispVelIm;
  const PetscPointJac hf1_veldisp = NULL;
  const PetscPointJac hf2_veldisp = NULL;
  const PetscPointJac hf3_veldisp = NULL;

  const PetscPointJac hf0_velvel = NULL;
  const PetscPointJac hf1_velvel = NULL;
  const PetscPointJac hf2_velvel = NULL;
  const PetscPointJac hf3_velvel = NULL;
    
  err = PetscDSSetJacobian(prob, i_disp, i_disp, hf0_dispdisp, hf1_dispdisp, hf2_dispdisp, hf3_dispdisp);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_disp, i_vel,  hf0_dispvel,  hf1_dispvel,  hf2_dispvel,  hf3_dispvel);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_disp, hf0_veldisp,  hf1_veldisp,  hf2_veldisp,  hf3_veldisp);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_vel,  hf0_velvel,   hf1_velvel,   hf2_velvel,   hf3_velvel);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobianImplicit


// ----------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,u,\dot{u}).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernelsLHSJacobianExplicit(const topology::Field& solution) const
{ // _setFEKernelsLHSJacobianExplicit
  PYLITH_METHOD_BEGIN;

  const PetscInt i_disp = solution.subfieldInfo("displacement").index;
  const PetscInt i_vel = solution.subfieldInfo("velocity").index;

  const PetscDM dm = solution.dmMesh();assert(dm);
  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(dm, &prob);PYLITH_CHECK_ERROR(err);

  // Jacobian kernels
  const PetscPointJac hf0_dispdisp = NULL;
  const PetscPointJac hf1_dispdisp = NULL;
  const PetscPointJac hf2_dispdisp = NULL;
  const PetscPointJac hf3_dispdisp = NULL;
  
  const PetscPointJac hf0_dispvel = (_useInertia) ? pylith_fekernels_hf0_dispdisp_ElasticityInertiaEx : NULL;
  const PetscPointJac hf1_dispvel = NULL;
  const PetscPointJac hf2_dispvel = NULL;
  const PetscPointJac hf3_dispvel = NULL;

  const PetscPointJac hf0_veldisp = pylith_fekernels_hf0_veldisp_DispVelEx;
  const PetscPointJac hf1_veldisp = NULL;
  const PetscPointJac hf2_veldisp = NULL;
  const PetscPointJac hf3_veldisp = NULL;

  const PetscPointJac hf0_velvel = NULL;
  const PetscPointJac hf1_velvel = NULL;
  const PetscPointJac hf2_velvel = NULL;
  const PetscPointJac hf3_velvel = NULL;
    
  err = PetscDSSetJacobian(prob, i_disp, i_disp, hf0_dispdisp, hf1_dispdisp, hf2_dispdisp, hf3_dispdisp);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_disp, i_vel,  hf0_dispvel,  hf1_dispvel,  hf2_dispvel,  hf3_dispvel);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_disp, hf0_veldisp,  hf1_veldisp,  hf2_veldisp,  hf3_veldisp);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_vel,  hf0_velvel,   hf1_velvel,   hf2_velvel,   hf3_velvel);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobianExplicit




// End of file 
