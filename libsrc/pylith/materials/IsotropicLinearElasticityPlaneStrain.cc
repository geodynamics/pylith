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
#include "pylith/fekernels/dispvel.h" // USES DispVel kernels
#include "pylith/fekernels/elasticity.h" // USES Elasticity kernels
#include "pylith/fekernels/linearelasticityplanestrain.h" // USES IsotropicLinearElasticityPlaneStrain kernels
}

#include "petscds.h"

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearElasticityPlaneStrain::IsotropicLinearElasticityPlaneStrain(void) :
  MaterialNew(2),
  _useInertia(false),
  _useBodyForce(false),
  _useInitialState(false)
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
// Use initial (reference) stress and strain in computation of stress
// and strain?
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::useInitialState(const bool value)
{ // useInitialState
  _useInitialState = value;
} // useInitialState


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

  // :ATTENTION: The order for subfieldAdd() must match the order of the auxiliary fields in the FE kernels.
  // Field 0: density
  const char* densityComponents[1] = {"density"};
  _auxFields->subfieldAdd("density", densityComponents, 1, pylith::topology::Field::SCALAR, this->auxFieldDiscretization("density"), densityScale, pylith::topology::FieldQuery::validatorPositive);
  _auxFieldsQuery->queryFn("density", pylith::topology::FieldQuery::dbQueryGeneric);

  // Field 1: shearModulus
  const char* shearModulusComponents[1] = {"shear_modulus"};
  _auxFields->subfieldAdd("shear_modulus", shearModulusComponents, 1, topology::Field::SCALAR, this->auxFieldDiscretization("shear_modulus"), pressureScale);
  _auxFieldsQuery->queryFn("shear_modulus", pylith::materials::Query::dbQueryShearModulus2D);

  // Field 2: bulkModulus
  const char* bulkModulusComponents[1] = {"bulk_modulus"};
  _auxFields->subfieldAdd("bulk_modulus", bulkModulusComponents, 1, topology::Field::SCALAR, this->auxFieldDiscretization("bulk_modulus"), pressureScale);
  _auxFieldsQuery->queryFn("bulk_modulus", pylith::materials::Query::dbQueryBulkModulus2D);

  // Field 3: body force
  if (_useBodyForce) {
    assert(2 == dimension());
    const char* components[2] = {"body_force_x", "body_force_y"};
    _auxFields->subfieldAdd("body_force", components, dimension(), topology::Field::VECTOR, this->auxFieldDiscretization("body_force"), forceScale);
    _auxFieldsQuery->queryFn("body_force", pylith::topology::FieldQuery::dbQueryGeneric);
  } // if

  // Fields 4 and 5: initial stress and initial strain
  if (_useInitialState) {
    const PylithInt stressSize = 4;
    const char* componentsStress[stressSize] = {"stress_xx", "stress_yy", "stress_xy", "stress_zz"};
    _auxFields->subfieldAdd("initial_stress", componentsStress, stressSize, topology::Field::OTHER, this->auxFieldDiscretization("initial_stress"), pressureScale);
    _auxFieldsQuery->queryFn("initial_stress", pylith::topology::FieldQuery::dbQueryGeneric);

    const PylithInt strainSize = 4;
    const char* componentsStrain[strainSize] = {"strain_xx", "strain_yy", "strain_xy", "strain_zz"};
    _auxFields->subfieldAdd("initial_strain", componentsStress, stressSize, topology::Field::OTHER, this->auxFieldDiscretization("initial_strain"), 1.0);
    _auxFieldsQuery->queryFn("initial_strain", pylith::topology::FieldQuery::dbQueryGeneric);
  } // if

  PYLITH_METHOD_END;
} // _auxFieldsSetup

// ----------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernelsRHSResidual(const topology::Field& solution) const
{ // _setFEKernelsRHSResidual
  PYLITH_METHOD_BEGIN;

  const PetscInt i_disp = solution.subfieldInfo("displacement").index;
  const PetscInt i_vel = solution.subfieldInfo("velocity").index;

  // Displacement
  const PetscPointFunc g0_u = (_useBodyForce) ? pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g0 : NULL;
  const PetscPointFunc g1_u = (!_useInitialState) ? pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1 : pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1_initstate;

  // Velocity
  const PetscPointFunc g0_v = pylith_fekernels_DispVel_g0;
  const PetscPointFunc g1_v = NULL;

  const PetscDM dm = solution.dmMesh();assert(dm);
  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(dm, &prob);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetResidual(prob, i_disp, g0_u, g1_u);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetResidual(prob, i_vel,  g0_v, g1_v);PYLITH_CHECK_ERROR(err);

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

  // Jacobian kernels
  const PetscPointJac Jg0_uu = NULL;
  const PetscPointJac Jg1_uu = NULL;
  const PetscPointJac Jg2_uu = NULL;
  const PetscPointJac Jg3_uu = pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jg3_uu;
  
  const PetscPointJac Jg0_uv = NULL;
  const PetscPointJac Jg1_uv = NULL;
  const PetscPointJac Jg2_uv = NULL;
  const PetscPointJac Jg3_uv = NULL;

  const PetscPointJac Jg0_vu = NULL;
  const PetscPointJac Jg1_vu = NULL;
  const PetscPointJac Jg2_vu = NULL;
  const PetscPointJac Jg3_vu = NULL;

  const PetscPointJac Jg0_vv = pylith_fekernels_DispVel_Jg0_vv;
  const PetscPointJac Jg1_vv = NULL;
  const PetscPointJac Jg2_vv = NULL;
  const PetscPointJac Jg3_vv = NULL;
    
  const PetscDM dm = solution.dmMesh();assert(dm);
  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(dm, &prob);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_disp, i_disp, Jg0_uu, Jg1_uu, Jg2_uu, Jg3_uu);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jg0_uv, Jg1_uv, Jg2_uv, Jg3_uv);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jg0_vu, Jg1_vu, Jg2_vu, Jg3_vu);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jg0_vv, Jg1_vv, Jg2_vv, Jg3_vv);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobian


// ----------------------------------------------------------------------
// Set kernels for LHS residual F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernelsLHSResidual(const topology::Field& solution) const
{ // _setFEKernelsLHSResidual
  PYLITH_METHOD_BEGIN;

  const PetscInt i_disp = solution.subfieldInfo("displacement").index;
  const PetscInt i_vel = solution.subfieldInfo("velocity").index;

  // Displacement
  const PetscPointFunc f0_u = (_useInertia) ? pylith_fekernels_IsotropicLinearElasticityPlaneStrain_f0 : NULL;
  const PetscPointFunc f1_u = NULL;

  // Velocity
  const PetscPointFunc f0_v = pylith_fekernels_DispVel_f0;
  const PetscPointFunc f1_v = NULL;

  const PetscDM dm = solution.dmMesh();assert(dm);
  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(dm, &prob);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetResidual(prob, i_disp, f0_u, f1_u);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetResidual(prob, i_vel,  f0_v, f1_v);PYLITH_CHECK_ERROR(err);

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

  // Jacobian kernels
  const PetscPointJac Jf0_uu = NULL;
  const PetscPointJac Jf1_uu = NULL;
  const PetscPointJac Jf2_uu = NULL;
  const PetscPointJac Jf3_uu = NULL;
  
  const PetscPointJac Jf0_uv = (_useInertia) ? pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0_uv_implicit : NULL;
  const PetscPointJac Jf1_uv = NULL;
  const PetscPointJac Jf2_uv = NULL;
  const PetscPointJac Jf3_uv = NULL;

  const PetscPointJac Jf0_vu = pylith_fekernels_DispVel_Jf0_vu_implicit;
  const PetscPointJac Jf1_vu = NULL;
  const PetscPointJac Jf2_vu = NULL;
  const PetscPointJac Jf3_vu = NULL;

  const PetscPointJac Jf0_vv = NULL;
  const PetscPointJac Jf1_vv = NULL;
  const PetscPointJac Jf2_vv = NULL;
  const PetscPointJac Jf3_vv = NULL;
    
  const PetscDM dm = solution.dmMesh();assert(dm);
  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(dm, &prob);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_disp, i_disp, Jf0_uu, Jf1_uu, Jf2_uu, Jf3_uu);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jf0_uv, Jf1_uv, Jf2_uv, Jf3_uv);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jf0_vu, Jf1_vu, Jf2_vu, Jf3_vu);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jf0_vv, Jf1_vv, Jf2_vv, Jf3_vv);PYLITH_CHECK_ERROR(err);

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

  // Jacobian kernels
  const PetscPointJac Jf0_uu = NULL;
  const PetscPointJac Jf1_uu = NULL;
  const PetscPointJac Jf2_uu = NULL;
  const PetscPointJac Jf3_uu = NULL;
  
  const PetscPointJac Jf0_uv = (_useInertia) ? pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0_uv_explicit : NULL;
  const PetscPointJac Jf1_uv = NULL;
  const PetscPointJac Jf2_uv = NULL;
  const PetscPointJac Jf3_uv = NULL;

  const PetscPointJac Jf0_vu = pylith_fekernels_DispVel_Jf0_vu_explicit;
  const PetscPointJac Jf1_vu = NULL;
  const PetscPointJac Jf2_vu = NULL;
  const PetscPointJac Jf3_vu = NULL;

  const PetscPointJac Jf0_vv = NULL;
  const PetscPointJac Jf1_vv = NULL;
  const PetscPointJac Jf2_vv = NULL;
  const PetscPointJac Jf3_vv = NULL;
    
  const PetscDM dm = solution.dmMesh();assert(dm);
  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(dm, &prob);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_disp, i_disp, Jf0_uu, Jf1_uu, Jf2_uu, Jf3_uu);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jf0_uv, Jf1_uv, Jf2_uv, Jf3_uv);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jf0_vu, Jf1_vu, Jf2_vu, Jf3_vu);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jf0_vv, Jf1_vv, Jf2_vv, Jf3_vv);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobianExplicit


// End of file 
