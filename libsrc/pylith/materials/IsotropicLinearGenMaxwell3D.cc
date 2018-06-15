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

#include "pylith/materials/IsotropicLinearGenMaxwell3D.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery

#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/Elasticity3D.hh" // USES Elasticity3D kernels
#include "pylith/fekernels/IsotropicLinearGenMaxwell3D.hh" // USES IsotropicLinearGenMaxwell3D kernels
#include "pylith/fekernels/DispVel.hh" // USES DispVel kernels

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petscds.h"

// ----------------------------------------------------------------------
const char* pylith::materials::IsotropicLinearGenMaxwell3D::_pyreComponent = "isotropiclineargenmaxwell3d";

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearGenMaxwell3D::IsotropicLinearGenMaxwell3D(void) :
    pylith::materials::Material(3),
    _useInertia(false),
    _useBodyForce(false),
    _useReferenceState(false)
{ // constructor
    pylith::utils::PyreComponent::name(_pyreComponent);
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearGenMaxwell3D::~IsotropicLinearGenMaxwell3D(void) {}


// ----------------------------------------------------------------------
// Include inertia?
void
pylith::materials::IsotropicLinearGenMaxwell3D::useInertia(const bool value) {
    PYLITH_COMPONENT_DEBUG("useInertia(value="<<value<<")");

    _useInertia = value;
} // useInertia


// ----------------------------------------------------------------------
// Include inertia?
bool
pylith::materials::IsotropicLinearGenMaxwell3D::useInertia(void) const {
    return _useInertia;
} // useInertia


// ----------------------------------------------------------------------
// Include body force?
void
pylith::materials::IsotropicLinearGenMaxwell3D::useBodyForce(const bool value) {
    PYLITH_COMPONENT_DEBUG("useBodyForce(value="<<value<<")");

    _useBodyForce = value;
} // useBodyForce


// ----------------------------------------------------------------------
// Include body force?
bool
pylith::materials::IsotropicLinearGenMaxwell3D::useBodyForce(void) const {
    return _useBodyForce;
} // useBodyForce


// ----------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void
pylith::materials::IsotropicLinearGenMaxwell3D::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ----------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool
pylith::materials::IsotropicLinearGenMaxwell3D::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::materials::IsotropicLinearGenMaxwell3D::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    // Verify solution contains expected fields.
    if (!solution.hasSubfield("displacement")) {
        throw std::runtime_error("Cannot find 'displacement' field in solution; required for material 'IsotropicLinearGenMaxwell3D'.");
    } // if
    if (_useInertia && !solution.hasSubfield("velocity")) {
        throw std::runtime_error("Cannot find 'velocity' field in solution; required for material 'IsotropicLinearGenMaxwell3D' with inertia.");
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Preinitialize material. Set names/sizes of auxiliary subfields.
void
pylith::materials::IsotropicLinearGenMaxwell3D::_auxFieldSetup(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxFieldSetup()");

    const int dim = 3;

    assert(_auxMaterialFactory);
    assert(_normalizer);
    _auxMaterialFactory->initialize(_auxField, *_normalizer, dim);

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.
    // NOTE:  At present, we assume exactly 3 Maxwell elements.

    _auxMaterialFactory->density(); // 0
    _auxMaterialFactory->shearModulus(); // 1
    _auxMaterialFactory->bulkModulus(); // 2
    _auxMaterialFactory->maxwellTimeGeneralizedMaxwell(); // 3
    _auxMaterialFactory->shearModulusRatioGeneralizedMaxwell(); // 4
    _auxMaterialFactory->viscousStrainGeneralizedMaxwell(); // 5
    _auxMaterialFactory->totalStrain(); // 6
    if (_gravityField) {
        _auxMaterialFactory->gravityField(_gravityField);
    } // if
    if (_useBodyForce) {
        _auxMaterialFactory->bodyForce();
    } // if
    if (_useReferenceState) {
        _auxMaterialFactory->referenceStress(); // numA-2
        _auxMaterialFactory->referenceStrain(); // numA-1
    } // if

    _updateStateVarsKernels["viscous_strain"] = pylith::fekernels::IsotropicLinearGenMaxwell3D::updateViscousStrain;
    _updateStateVarsKernels["total_strain"] = pylith::fekernels::IsotropicLinearGenMaxwell3D::updateTotalStrain;

    PYLITH_METHOD_END;
} // _auxFieldSetup

// ----------------------------------------------------------------------
// Set constants for problem.
void
pylith::materials::IsotropicLinearGenMaxwell3D::_setFEConstants(const pylith::topology::Field& solution,
																const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEConstants(solution="<<solution.label()<<", dt="<<dt<<")");

    const PetscInt numConstants = 1;
    PetscScalar constants[1];
    constants[0] = dt;

    const PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err); assert(prob);
    err = PetscDSSetConstants(prob, numConstants, constants); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setConstants

// ----------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::materials::IsotropicLinearGenMaxwell3D::_setFEKernelsRHSResidual(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsRHSResidual(solution="<<solution.label()<<")");

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;

    if (!solution.hasSubfield("velocity")) {
        // Displacement
        const PetscPointFunc g0u = (_gravityField && _useBodyForce) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::g0v_gravbodyforce :
                                   (_gravityField) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::g0v_grav :
                                   (_useBodyForce) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::g0v_bodyforce :
                                   NULL;
        const PetscPointFunc g1u = (!_useReferenceState) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::g1v : pylith::fekernels::IsotropicLinearGenMaxwell3D::g1v_refstate;

        err = PetscDSSetResidual(prob, i_disp, g0u, g1u); PYLITH_CHECK_ERROR(err);
    } else {
        const PetscInt i_vel = solution.subfieldInfo("velocity").index;

        // Displacement
        const PetscPointFunc g0u = pylith::fekernels::DispVel::g0u;
        const PetscPointFunc g1u = NULL;

        // Velocity
        const PetscPointFunc g0v = (_gravityField && _useBodyForce) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::g0v_gravbodyforce :
                                   (_gravityField) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::g0v_grav :
                                   (_useBodyForce) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::g0v_bodyforce :
                                   NULL;
        const PetscPointFunc g1v = (!_useReferenceState) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::g1v : pylith::fekernels::IsotropicLinearGenMaxwell3D::g1v_refstate;

        err = PetscDSSetResidual(prob, i_disp, g0u, g1u); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetResidual(prob, i_vel,  g0v, g1v); PYLITH_CHECK_ERROR(err);
    } // if/else


    PYLITH_METHOD_END;
} // _setFEKernelsRHSResidual


// ----------------------------------------------------------------------
// Set kernels for RHS Jacobian G(t,s).
void
pylith::materials::IsotropicLinearGenMaxwell3D::_setFEKernelsRHSJacobian(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsRHSJacobian(solution="<<solution.label()<<")");

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;

    if (!solution.hasSubfield("velocity")) {
        const PetscPointJac Jg0uu = NULL;
        const PetscPointJac Jg1uu = NULL;
        const PetscPointJac Jg2uu = NULL;
        const PetscPointJac Jg3uu = pylith::fekernels::IsotropicLinearGenMaxwell3D::Jg3vu;

        err = PetscDSSetJacobian(prob, i_disp, i_disp, Jg0uu, Jg1uu, Jg2uu, Jg3uu); PYLITH_CHECK_ERROR(err);
    } else {
        const PetscInt i_vel = solution.subfieldInfo("velocity").index;

        // Jacobian kernels
        const PetscPointJac Jg0uu = NULL;
        const PetscPointJac Jg1uu = NULL;
        const PetscPointJac Jg2uu = NULL;
        const PetscPointJac Jg3uu = NULL;

        const PetscPointJac Jg0uv = pylith::fekernels::DispVel::Jg0uv;
        const PetscPointJac Jg1uv = NULL;
        const PetscPointJac Jg2uv = NULL;
        const PetscPointJac Jg3uv = NULL;

        const PetscPointJac Jg0vu = NULL;
        const PetscPointJac Jg1vu = NULL;
        const PetscPointJac Jg2vu = NULL;
        const PetscPointJac Jg3vu = pylith::fekernels::IsotropicLinearGenMaxwell3D::Jg3vu;

        const PetscPointJac Jg0vv = NULL;
        const PetscPointJac Jg1vv = NULL;
        const PetscPointJac Jg2vv = NULL;
        const PetscPointJac Jg3vv = NULL;

        err = PetscDSSetJacobian(prob, i_disp, i_disp, Jg0uu, Jg1uu, Jg2uu, Jg3uu); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jg0uv, Jg1uv, Jg2uv, Jg3uv); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jg0vu, Jg1vu, Jg2vu, Jg3vu); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jg0vv, Jg1vv, Jg2vv, Jg3vv); PYLITH_CHECK_ERROR(err);

    } // if/else

    PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobian


// ----------------------------------------------------------------------
// Set kernels for LHS residual F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearGenMaxwell3D::_setFEKernelsLHSResidual(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsLHSResidual(solution="<<solution.label()<<")");

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;

    if (!solution.hasSubfield("velocity")) {
        // F(t,s,\dot{s}) = \vec{0}.
    } else {
        const PetscInt i_vel = solution.subfieldInfo("velocity").index;

        // Displacement
        const PetscPointFunc f0u = pylith::fekernels::DispVel::f0u;
        const PetscPointFunc f1u = NULL;

        // Velocity
        const PetscPointFunc f0v = (_useInertia) ? pylith::fekernels::DispVel::f0v : NULL;
        const PetscPointFunc f1v = NULL;

        err = PetscDSSetResidual(prob, i_disp, f0u, f1u); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetResidual(prob, i_vel,  f0v, f1v); PYLITH_CHECK_ERROR(err);
    } // if/else

    PYLITH_METHOD_END;
} // _setFEKernelsLHSResidual


// ----------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearGenMaxwell3D::_setFEKernelsLHSJacobian(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsLHSJacobian(solution="<<solution.label()<<")");

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;

    if (!solution.hasSubfield("velocity")) {
        // Jacobian kernels
        const PetscPointJac Jf0uu = pylith::fekernels::DispVel::Jf0uu_zero;
        const PetscPointJac Jf1uu = NULL;
        const PetscPointJac Jf2uu = NULL;
        const PetscPointJac Jf3uu = NULL;

        err = PetscDSSetJacobian(prob, i_disp, i_disp, Jf0uu, Jf1uu, Jf2uu, Jf3uu); PYLITH_CHECK_ERROR(err);
    } else {
        const PetscInt i_vel = solution.subfieldInfo("velocity").index;

        // Jacobian kernels
        const PetscPointJac Jf0uu = pylith::fekernels::DispVel::Jf0uu_utshift;
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

        const PetscPointJac Jf0vv = (_useInertia) ? pylith::fekernels::DispVel::Jf0uu_utshift : NULL;
        const PetscPointJac Jf1vv = NULL;
        const PetscPointJac Jf2vv = NULL;
        const PetscPointJac Jf3vv = NULL;

        err = PetscDSSetJacobian(prob, i_disp, i_disp, Jf0uu, Jf1uu, Jf2uu, Jf3uu); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jf0uv, Jf1uv, Jf2uv, Jf3uv); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jf0vu, Jf1vu, Jf2vu, Jf3vu); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jf0vv, Jf1vv, Jf2vv, Jf3vv); PYLITH_CHECK_ERROR(err);
    } // if/else

    PYLITH_METHOD_END;
} // _setFEKernelsLHSJacobian


// End of file
