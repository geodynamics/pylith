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

#include "IsotropicLinearIncompElasticityPlaneStrain.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/materials/Query.hh" // USES Query
extern "C" {
#include "pylith/fekernels/dispvel.h" // USES DispVel kernels
#include "pylith/fekernels/pressure.h" // USES Pressure kernels
#include "pylith/fekernels/elasticity.h" // USES Elasticity kernels
#include "pylith/fekernels/linearelasticityplanestrain.h" // USES IsotropicLinearElasticityPlaneStrain kernels
}

#include "petscds.h"

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearIncompElasticityPlaneStrain::IsotropicLinearIncompElasticityPlaneStrain(void) :
    MaterialNew(2),
    _useInertia(false),
    _useBodyForce(false),
    _useInitialState(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearIncompElasticityPlaneStrain::~IsotropicLinearIncompElasticityPlaneStrain(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Include inertia?
void
pylith::materials::IsotropicLinearIncompElasticityPlaneStrain::useInertia(const bool value)
{ // useInertia
    _useInertia = value;
} // useInertia

// ----------------------------------------------------------------------
// Include body force?
void
pylith::materials::IsotropicLinearIncompElasticityPlaneStrain::useBodyForce(const bool value)
{ // useBodyForce
    _useBodyForce = value;
} // useBodyForce

// ----------------------------------------------------------------------
// Include initial stress/strain?
void
pylith::materials::IsotropicLinearIncompElasticityPlaneStrain::useInitialState(const bool value)
{ // useInitialState
    _useInitialState = value;
} // useInitialState

// ----------------------------------------------------------------------
// Preinitialize material. Set names/sizes of auxiliary fields.
void
pylith::materials::IsotropicLinearIncompElasticityPlaneStrain::_auxFieldsSetup(void)
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
    const pylith::topology::Field::DiscretizeInfo& densityFEInfo = this->auxFieldDiscretization("density");
    _auxFields->subfieldAdd("density", densityComponents, 1, pylith::topology::Field::SCALAR, densityFEInfo.basisOrder, densityFEInfo.quadOrder, densityFEInfo.isBasisContinuous, densityScale, pylith::topology::FieldQuery::validatorPositive);
    _auxFieldsQuery->queryFn("density", pylith::topology::FieldQuery::dbQueryGeneric);

    // Field 1: shearModulus
    const char* shearModulusComponents[1] = {"shear_modulus"};
    const pylith::topology::Field::DiscretizeInfo& shearModulusFEInfo = this->auxFieldDiscretization("shear_modulus");
    _auxFields->subfieldAdd("shear_modulus", shearModulusComponents, 1, topology::Field::SCALAR, shearModulusFEInfo.basisOrder, shearModulusFEInfo.quadOrder, shearModulusFEInfo.isBasisContinuous, pressureScale);
    _auxFieldsQuery->queryFn("shear_modulus", pylith::materials::Query::dbQueryShearModulus2D);

    // Field 2: bulkModulus
    const char* bulkModulusComponents[1] = {"bulk_modulus"};
    const pylith::topology::Field::DiscretizeInfo& bulkModulusFEInfo = this->auxFieldDiscretization("bulk_modulus");
    _auxFields->subfieldAdd("bulk_modulus", bulkModulusComponents, 1, topology::Field::SCALAR, bulkModulusFEInfo.basisOrder, bulkModulusFEInfo.quadOrder, bulkModulusFEInfo.isBasisContinuous, pressureScale);
    _auxFieldsQuery->queryFn("bulk_modulus", pylith::materials::Query::dbQueryBulkModulus2D);

    // Field 3: body force
    if (_useBodyForce) {
        const char* bodyForceComponents[2] = {"body_force_x", "body_force_y"};
        const pylith::topology::Field::DiscretizeInfo& bodyForceFEInfo = this->auxFieldDiscretization("body_force");
        _auxFields->subfieldAdd("body force", bodyForceComponents, dimension(), topology::Field::VECTOR, bodyForceFEInfo.basisOrder, bodyForceFEInfo.quadOrder, bodyForceFEInfo.isBasisContinuous, forceScale);
        _auxFieldsQuery->queryFn("body force", pylith::topology::FieldQuery::dbQueryGeneric);
    } // if

    // Fields 4 and 5: initial stress and strain
    if (_useInitialState) {
        const PylithInt stressSize = 4;
        const char* componentsStress[stressSize] = {"stress_xx", "stress_yy", "stress_xy", "stress_zz"};
        const pylith::topology::Field::DiscretizeInfo& stressFEInfo = this->auxFieldDiscretization("reference_stress");
        _auxFields->subfieldAdd("initial_stress", componentsStress, stressSize, topology::Field::OTHER, stressFEInfo.basisOrder, stressFEInfo.quadOrder, stressFEInfo.isBasisContinuous, pressureScale);
        _auxFieldsQuery->queryFn("initial_stress", pylith::topology::FieldQuery::dbQueryGeneric);

        const PylithInt strainSize = 4;
        const char* componentsStrain[strainSize] = {"strain_xx", "strain_yy", "strain_xy", "strain_zz"};
        const pylith::topology::Field::DiscretizeInfo& strainFEInfo = this->auxFieldDiscretization("reference_strain");
        _auxFields->subfieldAdd("initial_strain", componentsStrain, strainSize, topology::Field::OTHER, strainFEInfo.basisOrder, strainFEInfo.quadOrder, strainFEInfo.isBasisContinuous, 1.0);
        _auxFieldsQuery->queryFn("initial_strain", pylith::topology::FieldQuery::dbQueryGeneric);
    } // if

    PYLITH_METHOD_END;
} // _auxFieldsSetup

// ----------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::materials::IsotropicLinearIncompElasticityPlaneStrain::_setFEKernelsRHSResidual(const topology::Field& solution) const
{ // _setFEKernelsRHSResidual
    PYLITH_METHOD_BEGIN;

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;
    const PetscInt i_vel = solution.subfieldInfo("velocity").index;
    const PetscInt i_pres = solution.subfieldInfo("pressure").index;

    // Displacement
    const PetscPointFunc g0u = pylith_fekernels_DispVel_g0u;
    const PetscPointFunc g1u = NULL;

    // Velocity
    const PetscPointFunc g0v = (_useBodyForce) ? pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g0v : NULL;
    const PetscPointFunc g1v = (_useInitialState) ? pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g1v : pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1v_initstate;

    // Pressure
    const PetscPointFunc g0p = pylith_fekernels_Pressure_g0p;
    const PetscPointFunc g1p = NULL;

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetResidual(prob, i_disp, g0u, g1u); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetResidual(prob, i_vel,  g0v, g1v); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetResidual(prob, i_pres, g0p, g1p); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSResidual


// ----------------------------------------------------------------------
// Set kernels for RHS Jacobian G(t,s).
void
pylith::materials::IsotropicLinearIncompElasticityPlaneStrain::_setFEKernelsRHSJacobian(const topology::Field& solution) const
{ // _setFEKernelsRHSJacobian
    PYLITH_METHOD_BEGIN;

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;
    const PetscInt i_vel = solution.subfieldInfo("velocity").index;
    const PetscInt i_pres = solution.subfieldInfo("pressure").index;

    // Jacobian kernels
    const PetscPointJac Jg0uu = NULL;
    const PetscPointJac Jg1uu = NULL;
    const PetscPointJac Jg2uu = NULL;
    const PetscPointJac Jg3uu = NULL;

    const PetscPointJac Jg0uv = pylith_fekernels_DispVel_Jg0uv;
    const PetscPointJac Jg1uv = NULL;
    const PetscPointJac Jg2uv = NULL;
    const PetscPointJac Jg3uv = NULL;

    const PetscPointJac Jg0up = NULL;
    const PetscPointJac Jg1up = NULL;
    const PetscPointJac Jg2up = NULL;
    const PetscPointJac Jg3up = NULL;

    const PetscPointJac Jg0vu = NULL;
    const PetscPointJac Jg1vu = NULL;
    const PetscPointJac Jg2vu = NULL;
    const PetscPointJac Jg3vu = pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_Jg3vu;

    const PetscPointJac Jg0vv = NULL;
    const PetscPointJac Jg1vv = NULL;
    const PetscPointJac Jg2vv = NULL;
    const PetscPointJac Jg3vv = NULL;

    const PetscPointJac Jg0vp = NULL;
    const PetscPointJac Jg1vp = NULL;
    const PetscPointJac Jg2vp = pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_Jg2vp;
    const PetscPointJac Jg3vp = NULL;

    const PetscPointJac Jg0pu = NULL;
    const PetscPointJac Jg1pu = NULL;
    const PetscPointJac Jg2pu = NULL;
    const PetscPointJac Jg3pu = NULL;

    const PetscPointJac Jg0pv = NULL;
    const PetscPointJac Jg1pv = NULL;
    const PetscPointJac Jg2pv = NULL;
    const PetscPointJac Jg3pv = NULL;

    const PetscPointJac Jg0pp = pylith_fekernels_Pressure_Jg0pp;
    const PetscPointJac Jg1pp = NULL;
    const PetscPointJac Jg2pp = NULL;
    const PetscPointJac Jg3pp = NULL;

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp, i_disp, Jg0uu, Jg1uu, Jg2uu, Jg3uu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp,  i_vel, Jg0uv, Jg1uv, Jg2uv, Jg3uv); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp, i_pres, Jg0up, Jg1up, Jg2up, Jg3up); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob,  i_vel, i_disp, Jg0vu, Jg1vu, Jg2vu, Jg3vu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob,  i_vel,  i_vel, Jg0vv, Jg1vv, Jg2vv, Jg3vv); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob,  i_vel, i_pres, Jg0vp, Jg1vp, Jg2vp, Jg3vp); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_pres, i_disp, Jg0pu, Jg1pu, Jg2pu, Jg3pu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_pres,  i_vel, Jg0pv, Jg1pv, Jg2pv, Jg3pv); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_pres, i_pres, Jg0pp, Jg1pp, Jg2pp, Jg3pp); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobian


// ----------------------------------------------------------------------
// Set kernels for LHS residual F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearIncompElasticityPlaneStrain::_setFEKernelsLHSResidual(const topology::Field& solution) const
{ // _setFEKernelsLHSResidual
    PYLITH_METHOD_BEGIN;

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;
    const PetscInt i_vel = solution.subfieldInfo("velocity").index;
    const PetscInt i_pres = solution.subfieldInfo("pressure").index;

    // Displacement
    const PetscPointFunc f0u = pylith_fekernels_DispVel_f0u;
    const PetscPointFunc f1u = NULL;

    // Velocity
    const PetscPointFunc f0v = (_useInertia) ? pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_f0v : NULL;
    const PetscPointFunc f1v = NULL;

    // Pressure
    const PetscPointFunc f0p = NULL;
    const PetscPointFunc f1p = NULL;

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetResidual(prob, i_disp, f0u, f1u); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetResidual(prob, i_vel,  f0v, f1v); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetResidual(prob, i_pres, f0p, f1p); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsLHSResidual


// ----------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearIncompElasticityPlaneStrain::_setFEKernelsLHSJacobianImplicit(const topology::Field& solution) const
{ // _setFEKernelsLHSJacobianImplicit
    PYLITH_METHOD_BEGIN;

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;
    const PetscInt i_vel = solution.subfieldInfo("velocity").index;
    const PetscInt i_pres = solution.subfieldInfo("pressure").index;

    // Jacobian kernels
    const PetscPointJac Jf0uu = pylith_fekernels_DispVel_Jf0uu_implicit;
    const PetscPointJac Jf1uu = NULL;
    const PetscPointJac Jf2uu = NULL;
    const PetscPointJac Jf3uu = NULL;

    const PetscPointJac Jf0uv = NULL;
    const PetscPointJac Jf1uv = NULL;
    const PetscPointJac Jf2uv = NULL;
    const PetscPointJac Jf3uv = NULL;

    const PetscPointJac Jf0up = NULL;
    const PetscPointJac Jf1up = NULL;
    const PetscPointJac Jf2up = NULL;
    const PetscPointJac Jf3up = NULL;

    const PetscPointJac Jf0pu = NULL;
    const PetscPointJac Jf1pu = NULL;
    const PetscPointJac Jf2pu = NULL;
    const PetscPointJac Jf3pu = NULL;

    const PetscPointJac Jf0vu = NULL;
    const PetscPointJac Jf1vu = NULL;
    const PetscPointJac Jf2vu = NULL;
    const PetscPointJac Jf3vu = NULL;

    const PetscPointJac Jf0vv = (_useInertia) ? pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0vv_implicit : NULL;
    const PetscPointJac Jf1vv = NULL;
    const PetscPointJac Jf2vv = NULL;
    const PetscPointJac Jf3vv = NULL;

    const PetscPointJac Jf0vp = NULL;
    const PetscPointJac Jf1vp = NULL;
    const PetscPointJac Jf2vp = NULL;
    const PetscPointJac Jf3vp = NULL;

    const PetscPointJac Jf0pv = NULL;
    const PetscPointJac Jf1pv = NULL;
    const PetscPointJac Jf2pv = NULL;
    const PetscPointJac Jf3pv = NULL;

    const PetscPointJac Jf0pp = NULL;
    const PetscPointJac Jf1pp = NULL;
    const PetscPointJac Jf2pp = NULL;
    const PetscPointJac Jf3pp = NULL;

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp, i_disp, Jf0uu, Jf1uu, Jf2uu, Jf3uu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jf0uv, Jf1uv, Jf2uv, Jf3uv); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp, i_pres, Jf0up, Jf1up, Jf2up, Jf3up); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jf0vu, Jf1vu, Jf2vu, Jf3vu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jf0vv, Jf1vv, Jf2vv, Jf3vv); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_vel,  i_pres, Jf0vp, Jf1vp, Jf2vp, Jf3vp); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_pres, i_disp, Jf0pu, Jf1pu, Jf2pu, Jf3pu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_pres, i_vel,  Jf0pv, Jf1pv, Jf2pv, Jf3pv); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_pres, i_pres, Jf0pp, Jf1pp, Jf2pp, Jf3pp); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsLHSJacobianImplicit


// ----------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearIncompElasticityPlaneStrain::_setFEKernelsLHSJacobianExplicit(const topology::Field& solution) const
{ // _setFEKernelsLHSJacobianExplicit
    PYLITH_METHOD_BEGIN;

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;
    const PetscInt i_vel = solution.subfieldInfo("velocity").index;
    const PetscInt i_pres = solution.subfieldInfo("pressure").index;

    // Jacobian kernels
    const PetscPointJac Jf0uu = pylith_fekernels_DispVel_Jf0uu_explicit;
    const PetscPointJac Jf1uu = NULL;
    const PetscPointJac Jf2uu = NULL;
    const PetscPointJac Jf3uu = NULL;

    const PetscPointJac Jf0uv = NULL;
    const PetscPointJac Jf1uv = NULL;
    const PetscPointJac Jf2uv = NULL;
    const PetscPointJac Jf3uv = NULL;

    const PetscPointJac Jf0up = NULL;
    const PetscPointJac Jf1up = NULL;
    const PetscPointJac Jf2up = NULL;
    const PetscPointJac Jf3up = NULL;

    const PetscPointJac Jf0vu = NULL;
    const PetscPointJac Jf1vu = NULL;
    const PetscPointJac Jf2vu = NULL;
    const PetscPointJac Jf3vu = NULL;

    const PetscPointJac Jf0vv = (_useInertia) ? pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0vv_explicit : NULL;
    const PetscPointJac Jf1vv = NULL;
    const PetscPointJac Jf2vv = NULL;
    const PetscPointJac Jf3vv = NULL;

    const PetscPointJac Jf0vp = NULL;
    const PetscPointJac Jf1vp = NULL;
    const PetscPointJac Jf2vp = NULL;
    const PetscPointJac Jf3vp = NULL;

    const PetscPointJac Jf0pu = NULL;
    const PetscPointJac Jf1pu = NULL;
    const PetscPointJac Jf2pu = NULL;
    const PetscPointJac Jf3pu = NULL;

    const PetscPointJac Jf0pv = NULL;
    const PetscPointJac Jf1pv = NULL;
    const PetscPointJac Jf2pv = NULL;
    const PetscPointJac Jf3pv = NULL;

    const PetscPointJac Jf0pp = NULL;
    const PetscPointJac Jf1pp = NULL;
    const PetscPointJac Jf2pp = NULL;
    const PetscPointJac Jf3pp = NULL;

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp, i_disp, Jf0uu, Jf1uu, Jf2uu, Jf3uu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jf0uv, Jf1uv, Jf2uv, Jf3uv); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp, i_pres, Jf0up, Jf1up, Jf2up, Jf3up); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jf0vu, Jf1vu, Jf2vu, Jf3vu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jf0vv, Jf1vv, Jf2vv, Jf3vv); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_vel,  i_pres, Jf0vp, Jf1vp, Jf2vp, Jf3vp); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_pres, i_disp, Jf0pu, Jf1pu, Jf2pu, Jf3pu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_pres, i_vel,  Jf0pv, Jf1pv, Jf2pv, Jf3pv); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_pres, i_pres, Jf0pp, Jf1pp, Jf2pp, Jf3pp); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobianExplicit


// End of file
