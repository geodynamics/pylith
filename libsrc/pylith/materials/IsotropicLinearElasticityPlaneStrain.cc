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
// Set residual and Jacobian kernels.
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernels(const topology::Field& field) const
{ // _setFEKernels
  PYLITH_METHOD_BEGIN;

  const PetscInt i_disp = field.subfieldInfo("displacement").index;
  const PetscInt i_vel = (_useInertia) ? field.subfieldInfo("velocity").index : -1;

  PetscErrorCode err;

  const PetscDM dm = field.dmMesh();assert(dm);
  PetscDS prob = NULL;
  err = DMGetDS(dm, &prob);PYLITH_CHECK_ERROR(err);

  // Residual kernels
  const PetscPointFunc f0_disp = (_useInertia && _useBodyForce) ? pylith_fekernels_f0_ElasticityInertiaBodyForce : 
    (_useInertia) ? pylith_fekernels_f0_ElasticityInertia : 
    (_useBodyForce) ? pylith_fekernels_f0_ElasticityBodyForce : NULL;
  const PetscPointFunc f1_disp = pylith_fekernels_f1_IsotropicLinearElasticityPlaneStrain;
  err = PetscDSSetResidual(prob, i_disp, f0_disp, f1_disp);PYLITH_CHECK_ERROR(err);

  if (_useInertia) {
    const PetscPointFunc f0_vel = pylith_fekernels_f0_DispVel;
    const PetscPointFunc f1_vel = NULL;
    err = PetscDSSetResidual(prob, i_vel, f0_vel, f1_vel);PYLITH_CHECK_ERROR(err);
  } // if

  // Jacobian kernels
  const PetscPointJac g0_dispdisp = NULL;
  const PetscPointJac g1_dispdisp = NULL;
  const PetscPointJac g2_dispdisp = NULL;
  const PetscPointJac g3_dispdisp = pylith_fekernels_g3_uu_IsotropicLinearElasticityPlaneStrain;
  err = PetscDSSetJacobian(prob, i_disp, i_disp, g0_dispdisp, g1_dispdisp, g2_dispdisp, g3_dispdisp);PYLITH_CHECK_ERROR(err);
  
  if (_useInertia) {
    const PetscPointJac g0_dispvel = pylith_fekernels_g0_uv_ElasticityInertia;
    const PetscPointJac g1_dispvel = NULL;
    const PetscPointJac g2_dispvel = NULL;
    const PetscPointJac g3_dispvel = NULL;
    const PetscPointJac g0_veldisp = pylith_fekernels_g0_vu_DispVel;
    const PetscPointJac g1_veldisp = NULL;
    const PetscPointJac g2_veldisp = NULL;
    const PetscPointJac g3_veldisp = NULL;
    const PetscPointJac g0_velvel = pylith_fekernels_g0_vv_DispVel;
    const PetscPointJac g1_velvel = NULL;
    const PetscPointJac g2_velvel = NULL;
    const PetscPointJac g3_velvel = NULL;
    
    err = PetscDSSetJacobian(prob, i_disp, i_vel, g0_dispvel, g1_dispvel, g2_dispvel, g3_dispvel);PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_vel, i_disp, g0_veldisp, g1_veldisp, g2_veldisp, g3_veldisp);PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_vel, i_vel, g0_velvel, g1_velvel, g2_velvel, g3_velvel);PYLITH_CHECK_ERROR(err);
  } // if

  PYLITH_METHOD_END;
} // _setFEKernels


// End of file 
