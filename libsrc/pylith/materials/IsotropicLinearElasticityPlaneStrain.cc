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

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include "petscds.h"

// ----------------------------------------------------------------------
const char* pylith::materials::IsotropicLinearElasticityPlaneStrain::_pyreComponent = "isotropiclinearelasticityplanestrain";

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearElasticityPlaneStrain::IsotropicLinearElasticityPlaneStrain(void) :
    MaterialNew(2),
    _useInertia(false),
    _useBodyForce(false),
    _useReferenceState(false)
{ // constructor
    PyreComponent::name(_pyreComponent);
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
    PYLITH_JOURNAL_DEBUG("useInertia(value="<<value<<")");

    _useInertia = value;
} // useInertia


// ----------------------------------------------------------------------
// Include inertia?
bool
pylith::materials::IsotropicLinearElasticityPlaneStrain::useInertia(void) const
{ // useInertia
    return _useInertia;
} // useInertia


// ----------------------------------------------------------------------
// Include body force?
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::useBodyForce(const bool value)
{ // useBodyForce
    PYLITH_JOURNAL_DEBUG("useBodyForce(value="<<value<<")");

    _useBodyForce = value;
} // useBodyForce


// ----------------------------------------------------------------------
// Include body force?
bool
pylith::materials::IsotropicLinearElasticityPlaneStrain::useBodyForce(void) const
{ // useBodyForce
    return _useBodyForce;
} // useBodyForce


// ----------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::useReferenceState(const bool value)
{ // useReferenceState
    PYLITH_JOURNAL_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ----------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool
pylith::materials::IsotropicLinearElasticityPlaneStrain::useReferenceState(void) const
{ // useReferenceState
    return _useReferenceState;
} // useReferenceState


// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::verifyConfiguration(const pylith::topology::Field& solution) const
{ // verifyConfiguration
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    // Verify solution contains expected fields.
    if (!solution.hasSubfield("displacement")) {
        throw std::runtime_error("Cannot find 'displacement' field in solution; required for material 'IsotropicLinearElasticityPlaneStrain'.");
    } // if
    if (_useInertia && !solution.hasSubfield("velocity")) {
        throw std::runtime_error("Cannot find 'velocity' field in solution; required for material 'IsotropicLinearElasticityPlaneStrain' with inertia.");
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Preinitialize material. Set names/sizes of auxiliary fields.
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_auxFieldsSetup(void)
{ // _auxFieldsSetup
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_auxFieldsSetup()");

    PYLITH_JOURNAL_ERROR(":TODO: Add auxiliary field for gravitational acceleration vector");

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
        assert(2 == dimension());
        const char* components[2] = {"body_force_x", "body_force_y"};
        const pylith::topology::Field::DiscretizeInfo& bodyForceFEInfo = this->auxFieldDiscretization("body_force");
        _auxFields->subfieldAdd("body_force", components, dimension(), topology::Field::VECTOR, bodyForceFEInfo.basisOrder, bodyForceFEInfo.quadOrder, bodyForceFEInfo.isBasisContinuous, forceScale);
        _auxFieldsQuery->queryFn("body_force", pylith::topology::FieldQuery::dbQueryGeneric);
    } // if

    // Fields 4 and 5: reference stress and reference strain
    if (_useReferenceState) {
        const PylithInt stressSize = 4;
        const char* componentsStress[stressSize] = {"stress_xx", "stress_yy", "stress_xy", "stress_zz"};
        const pylith::topology::Field::DiscretizeInfo& stressFEInfo = this->auxFieldDiscretization("reference_stress");
        _auxFields->subfieldAdd("reference_stress", componentsStress, stressSize, topology::Field::OTHER, stressFEInfo.basisOrder, stressFEInfo.quadOrder, stressFEInfo.isBasisContinuous, pressureScale);
        _auxFieldsQuery->queryFn("reference_stress", pylith::topology::FieldQuery::dbQueryGeneric);

        const PylithInt strainSize = 4;
        const char* componentsStrain[strainSize] = {"strain_xx", "strain_yy", "strain_xy", "strain_zz"};
        const pylith::topology::Field::DiscretizeInfo& strainFEInfo = this->auxFieldDiscretization("reference_strain");
        _auxFields->subfieldAdd("reference_strain", componentsStrain, strainSize, topology::Field::OTHER, strainFEInfo.basisOrder, strainFEInfo.quadOrder, strainFEInfo.isBasisContinuous, 1.0);
        _auxFieldsQuery->queryFn("reference_strain", pylith::topology::FieldQuery::dbQueryGeneric);
    } // if

    PYLITH_METHOD_END;
} // _auxFieldsSetup

// ----------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernelsRHSResidual(const topology::Field& solution) const
{ // _setFEKernelsRHSResidual
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_setFEKernelsRHSResidual(solution="<<solution.label()<<")");

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;
    const PetscInt i_vel = solution.subfieldInfo("velocity").index;

    // Displacement
    const PetscPointFunc g0u = pylith_fekernels_DispVel_g0u;
    const PetscPointFunc g1u = NULL;

    // Velocity
    const PetscPointFunc g0v = (_useBodyForce) ? pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g0v : NULL;
    const PetscPointFunc g1v = (!_useReferenceState) ? pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1v : pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1v_initstate;

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetResidual(prob, i_disp, g0u, g1u); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetResidual(prob, i_vel,  g0v, g1v); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSResidual


// ----------------------------------------------------------------------
// Set kernels for RHS Jacobian G(t,s).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernelsRHSJacobian(const topology::Field& solution) const
{ // _setFEKernelsRHSJacobian
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_setFEKernelsRHSJacobian(solution="<<solution.label()<<")");

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;
    const PetscInt i_vel = solution.subfieldInfo("velocity").index;

    // Jacobian kernels
    const PetscPointJac Jg0uu = NULL;
    const PetscPointJac Jg1uu = NULL;
    const PetscPointJac Jg2uu = NULL;
    const PetscPointJac Jg3uu = NULL;

    const PetscPointJac Jg0uv = pylith_fekernels_DispVel_Jg0uv;
    const PetscPointJac Jg1uv = NULL;
    const PetscPointJac Jg2uv = NULL;
    const PetscPointJac Jg3uv = NULL;

    const PetscPointJac Jg0vu = NULL;
    const PetscPointJac Jg1vu = NULL;
    const PetscPointJac Jg2vu = NULL;
    const PetscPointJac Jg3vu = pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jg3vu;

    const PetscPointJac Jg0vv = NULL;
    const PetscPointJac Jg1vv = NULL;
    const PetscPointJac Jg2vv = NULL;
    const PetscPointJac Jg3vv = NULL;

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp, i_disp, Jg0uu, Jg1uu, Jg2uu, Jg3uu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jg0uv, Jg1uv, Jg2uv, Jg3uv); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jg0vu, Jg1vu, Jg2vu, Jg3vu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jg0vv, Jg1vv, Jg2vv, Jg3vv); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobian


// ----------------------------------------------------------------------
// Set kernels for LHS residual F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernelsLHSResidual(const topology::Field& solution) const
{ // _setFEKernelsLHSResidual
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_setFEKernelsLHSResidual(solution="<<solution.label()<<")");

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;
    const PetscInt i_vel = solution.subfieldInfo("velocity").index;

    // Displacement
    const PetscPointFunc f0u = pylith_fekernels_DispVel_f0u;
    const PetscPointFunc f1u = NULL;

    // Velocity
    const PetscPointFunc f0v = (_useInertia) ? pylith_fekernels_IsotropicLinearElasticityPlaneStrain_f0v : NULL;
    const PetscPointFunc f1v = NULL;

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetResidual(prob, i_disp, f0u, f1u); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetResidual(prob, i_vel,  f0v, f1v); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsLHSResidual


// ----------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernelsLHSJacobianImplicit(const topology::Field& solution) const
{ // _setFEKernelsLHSJacobianImplicit
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_setFEKernelsLHSJacobianImplicit(solution="<<solution.label()<<")");

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;
    const PetscInt i_vel = solution.subfieldInfo("velocity").index;

    // Jacobian kernels
    const PetscPointJac Jf0uu = pylith_fekernels_DispVel_Jf0uu_implicit;
    const PetscPointJac Jf1uu = NULL;
    const PetscPointJac Jf2uu = NULL;
    const PetscPointJac Jf3uu = NULL;

    const PetscPointJac Jf0uv = NULL;
    const PetscPointJac Jf1uv = NULL;
    const PetscPointJac Jf2uv = NULL;
    const PetscPointJac Jf3uv = NULL;

    const PetscPointJac Jf0vu = NULL;
    const PetscPointJac Jf1vu = NULL;
    const PetscPointJac Jf2vu = NULL;
    const PetscPointJac Jf3vu = NULL;

    const PetscPointJac Jf0vv = (_useInertia) ? pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0vv_implicit : NULL;
    const PetscPointJac Jf1vv = NULL;
    const PetscPointJac Jf2vv = NULL;
    const PetscPointJac Jf3vv = NULL;

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp, i_disp, Jf0uu, Jf1uu, Jf2uu, Jf3uu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jf0uv, Jf1uv, Jf2uv, Jf3uv); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jf0vu, Jf1vu, Jf2vu, Jf3vu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jf0vv, Jf1vv, Jf2vv, Jf3vv); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobianImplicit


// ----------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setFEKernelsLHSJacobianExplicit(const topology::Field& solution) const
{ // _setFEKernelsLHSJacobianExplicit
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_setFEKernelsLHSJacobianExplicit(solution="<<solution.label()<<")");

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;
    const PetscInt i_vel = solution.subfieldInfo("velocity").index;

    // Jacobian kernels
    const PetscPointJac Jf0uu = pylith_fekernels_DispVel_Jf0uu_explicit;
    const PetscPointJac Jf1uu = NULL;
    const PetscPointJac Jf2uu = NULL;
    const PetscPointJac Jf3uu = NULL;

    const PetscPointJac Jf0uv = NULL;
    const PetscPointJac Jf1uv = NULL;
    const PetscPointJac Jf2uv = NULL;
    const PetscPointJac Jf3uv = NULL;

    const PetscPointJac Jf0vu = NULL;
    const PetscPointJac Jf1vu = NULL;
    const PetscPointJac Jf2vu = NULL;
    const PetscPointJac Jf3vu = NULL;

    const PetscPointJac Jf0vv = (_useInertia) ? pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0vv_explicit : NULL;
    const PetscPointJac Jf1vv = NULL;
    const PetscPointJac Jf2vv = NULL;
    const PetscPointJac Jf3vv = NULL;

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp, i_disp, Jf0uu, Jf1uu, Jf2uu, Jf3uu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jf0uv, Jf1uv, Jf2uv, Jf3uv); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jf0vu, Jf1vu, Jf2vu, Jf3vu); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jf0vv, Jf1vv, Jf2vv, Jf3vv); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobianExplicit


// End of file
