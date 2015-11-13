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
{ // _auxFieldsSetup
  PYLITH_METHOD_BEGIN;

  // Set subfields in auxiliary fields.
  assert(_normalizer);
  const PylithReal densityScale = _normalizer->densityScale();
  const PylithReal pressureScale = _normalizer->pressureScale();
  const PylithReal lengthScale = _normalizer->lengthScale();
  const PylithReal timeScale = _normalizer->timeScale();
  const PylithReal forceScale = densityScale * lengthScale / (timeScale * timeScale);

  // :ATTENTION: The order here must match the order of the auxiliary fields in the FE kernels.
  _auxFields->subfieldAdd("density", 1, topology::Field::SCALAR, this->auxFieldDiscretization("density"), densityScale);
  _auxFields->subfieldAdd("mu", 1, topology::Field::SCALAR, this->auxFieldDiscretization("mu"), pressureScale);
  _auxFields->subfieldAdd("lambda", 1, topology::Field::SCALAR, this->auxFieldDiscretization("lambda"), pressureScale);
  if (_useBodyForce) {
    _auxFields->subfieldAdd("body force", dimension(), topology::Field::VECTOR, this->auxFieldDiscretization("body force"), forceScale);
  } // if

  // Order does not matter.
  _auxFieldsQuery->queryFn("density", pylith::materials::Query::dbQueryDensity2D);
  _auxFieldsQuery->queryFn("mu", pylith::materials::Query::dbQueryMu2D);
  _auxFieldsQuery->queryFn("lambda", pylith::materials::Query::dbQueryLambda2D);
  if (_useBodyForce) {
    _auxFieldsQuery->queryFn("body force", pylith::materials::Query::dbQueryBodyForce2D);
  } // if

  PYLITH_METHOD_END;
} // _auxFieldsSetup

// ----------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernelsRHSResidual(const topology::Field& solution) const
{ // _setFEKernelsRHSResidual
  PYLITH_METHOD_BEGIN;

  const PetscDM dm = solution.dmMesh();assert(dm);
  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(dm, &prob);PYLITH_CHECK_ERROR(err);

  // Displacement
  const PetscInt i_disp = solution.subfieldInfo("displacement").index;
  const PetscPointFunc g0_disp = (_useBodyForce) ? pylith_fekernels_g0_Elasticity : NULL;
  //const PetscPointFunc g1_disp = pylith_fekernels_g1_IsotropicLinearElasticityPlaneStrain;
  //err = PetscDSSetResidual(prob, i_disp, g0_disp, g1_disp);PYLITH_CHECK_ERROR(err);

  // Velocity
  const PetscInt i_vel = solution.subfieldInfo("velocity").index;
  const PetscPointFunc g0_vel = pylith_fekernels_g0_DispVel;
  const PetscPointFunc g1_vel = NULL;
  err = PetscDSSetResidual(prob, i_vel, g0_vel, g1_vel);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _setFEKernelsRHSResidual


// ----------------------------------------------------------------------
// Set kernels for RHS Jacobian G(t,s).
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
  const PetscPointJac Jg0_dispdisp = NULL;
  const PetscPointJac Jg1_dispdisp = NULL;
  const PetscPointJac Jg2_dispdisp = NULL;
  //const PetscPointJac Jg3_dispdisp = pylith_fekernels_Jg3_dispdisp_IsotropicLinearElasticityPlaneStrain;
  
  const PetscPointJac Jg0_dispvel = NULL;
  const PetscPointJac Jg1_dispvel = NULL;
  const PetscPointJac Jg2_dispvel = NULL;
  const PetscPointJac Jg3_dispvel = NULL;

  const PetscPointJac Jg0_veldisp = NULL;
  const PetscPointJac Jg1_veldisp = NULL;
  const PetscPointJac Jg2_veldisp = NULL;
  const PetscPointJac Jg3_veldisp = NULL;

  const PetscPointJac Jg0_velvel = pylith_fekernels_Jg0_velvel_DispVel;
  const PetscPointJac Jg1_velvel = NULL;
  const PetscPointJac Jg2_velvel = NULL;
  const PetscPointJac Jg3_velvel = NULL;
    
  //err = PetscDSSetJacobian(prob, i_disp, i_disp, Jg0_dispdisp, Jg1_dispdisp, Jg2_dispdisp, Jg3_dispdisp);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jg0_dispvel,  Jg1_dispvel,  Jg2_dispvel,  Jg3_dispvel);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jg0_veldisp,  Jg1_veldisp,  Jg2_veldisp,  Jg3_veldisp);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jg0_velvel,   Jg1_velvel,   Jg2_velvel,   Jg3_velvel);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobian


// ----------------------------------------------------------------------
// Set kernels for LHS residual F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernelsLHSResidual(const topology::Field& solution) const
{ // _setFEKernelsLHSResidual
  PYLITH_METHOD_BEGIN;

  const PetscDM dm = solution.dmMesh();assert(dm);
  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(dm, &prob);PYLITH_CHECK_ERROR(err);

  // Displacement
  const PetscInt i_disp = solution.subfieldInfo("displacement").index;
  const PetscPointFunc f0_disp = (_useInertia) ? pylith_fekernels_f0_Elasticity : NULL;
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
// Set kernels for LHS Jacobian F(t,s,\dot{s}).
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
  const PetscPointJac Jf0_dispdisp = NULL;
  const PetscPointJac Jf1_dispdisp = NULL;
  const PetscPointJac Jf2_dispdisp = NULL;
  const PetscPointJac Jf3_dispdisp = NULL;
  
  //const PetscPointJac Jf0_dispvel = (_useInertia) ? pylith_fekernels_Jf0_dispdisp_ElasticityImplicit : NULL;
  const PetscPointJac Jf1_dispvel = NULL;
  const PetscPointJac Jf2_dispvel = NULL;
  const PetscPointJac Jf3_dispvel = NULL;

  //const PetscPointJac Jf0_veldisp = pylith_fekernels_Jf0_veldisp_DispVelImplicit;
  const PetscPointJac Jf1_veldisp = NULL;
  const PetscPointJac Jf2_veldisp = NULL;
  const PetscPointJac Jf3_veldisp = NULL;

  const PetscPointJac Jf0_velvel = NULL;
  const PetscPointJac Jf1_velvel = NULL;
  const PetscPointJac Jf2_velvel = NULL;
  const PetscPointJac Jf3_velvel = NULL;
    
  err = PetscDSSetJacobian(prob, i_disp, i_disp, Jf0_dispdisp, Jf1_dispdisp, Jf2_dispdisp, Jf3_dispdisp);PYLITH_CHECK_ERROR(err);
  //err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jf0_dispvel,  Jf1_dispvel,  Jf2_dispvel,  Jf3_dispvel);PYLITH_CHECK_ERROR(err);
  //err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jf0_veldisp,  Jf1_veldisp,  Jf2_veldisp,  Jf3_veldisp);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jf0_velvel,   Jf1_velvel,   Jf2_velvel,   Jf3_velvel);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobianImplicit


// ----------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,s,\dot{s}).
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
  const PetscPointJac Jf0_dispdisp = NULL;
  const PetscPointJac Jf1_dispdisp = NULL;
  const PetscPointJac Jf2_dispdisp = NULL;
  const PetscPointJac Jf3_dispdisp = NULL;
  
  //const PetscPointJac Jf0_dispvel = (_useInertia) ? pylith_fekernels_Jf0_dispdisp_ElasticityExplicit : NULL;
  const PetscPointJac Jf1_dispvel = NULL;
  const PetscPointJac Jf2_dispvel = NULL;
  const PetscPointJac Jf3_dispvel = NULL;

  const PetscPointJac Jf0_veldisp = pylith_fekernels_Jf0_veldisp_DispVelExplicit;
  const PetscPointJac Jf1_veldisp = NULL;
  const PetscPointJac Jf2_veldisp = NULL;
  const PetscPointJac Jf3_veldisp = NULL;

  const PetscPointJac Jf0_velvel = NULL;
  const PetscPointJac Jf1_velvel = NULL;
  const PetscPointJac Jf2_velvel = NULL;
  const PetscPointJac Jf3_velvel = NULL;
    
  err = PetscDSSetJacobian(prob, i_disp, i_disp, Jf0_dispdisp, Jf1_dispdisp, Jf2_dispdisp, Jf3_dispdisp);PYLITH_CHECK_ERROR(err);
  //err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jf0_dispvel,  Jf1_dispvel,  Jf2_dispvel,  Jf3_dispvel);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jf0_veldisp,  Jf1_veldisp,  Jf2_veldisp,  Jf3_veldisp);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jf0_velvel,   Jf1_velvel,   Jf2_velvel,   Jf3_velvel);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobianExplicit




// End of file 
