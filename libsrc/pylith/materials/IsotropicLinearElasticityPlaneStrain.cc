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

#include "pylith/materials/IsotropicLinearElasticityPlaneStrain.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryElastic.hh" // USES AuxiliaryFactoryElastic
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/ElasticityPlaneStrain.hh" // USES ElasticityPlaneStrain kernels
#include "pylith/fekernels/IsotropicLinearElasticityPlaneStrain.hh" // USES IsotropicLinearElasticityPlaneStrain kernels
#include "pylith/fekernels/DispVel.hh" // USES DispVel kernels

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::Integrator::ResidualKernels ResidualKernels;
typedef pylith::feassemble::Integrator::JacobianKernels JacobianKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearElasticityPlaneStrain::IsotropicLinearElasticityPlaneStrain(void) :
    pylith::materials::Material(2),
    _useInertia(false),
    _useBodyForce(false),
    _useReferenceState(false),
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryElastic) {
    pylith::utils::PyreComponent::name("isotropiclinearelasticityplanestrain");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearElasticityPlaneStrain::~IsotropicLinearElasticityPlaneStrain(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::deallocate(void) {
    Material::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Include inertia?
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::useInertia(const bool value) {
    PYLITH_COMPONENT_DEBUG("useInertia(value="<<value<<")");

    _useInertia = value;
} // useInertia


// ---------------------------------------------------------------------------------------------------------------------
// Include inertia?
bool
pylith::materials::IsotropicLinearElasticityPlaneStrain::useInertia(void) const {
    return _useInertia;
} // useInertia


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::useBodyForce(const bool value) {
    PYLITH_COMPONENT_DEBUG("useBodyForce(value="<<value<<")");

    _useBodyForce = value;
} // useBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
bool
pylith::materials::IsotropicLinearElasticityPlaneStrain::useBodyForce(void) const {
    return _useBodyForce;
} // useBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool
pylith::materials::IsotropicLinearElasticityPlaneStrain::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::materials::IsotropicLinearElasticityPlaneStrain::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;

    pylith::feassemble::IntegratorDomain* integrator = new pylith::feassemble::IntegratorDomain(this);assert(integrator);
    integrator->setMaterialId(getMaterialId());

    _setKernelsRHSResidual(integrator, solution);
    _setKernelsRHSJacobian(integrator, solution);
    _setKernelsLHSResidual(integrator, solution);
    _setKernelsLHSJacobian(integrator, solution);
    // No state variables.
    // _setKernelsDerivedFields(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::materials::IsotropicLinearElasticityPlaneStrain::createAuxiliaryField(const pylith::topology::Field& solution,
                                                                              const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->label("IsotropicLinearElasticityPlaneStrain auxiliary field");

    const int dim = 2;

    assert(_auxiliaryFactory);
    assert(_normalizer);
    _auxiliaryFactory->initialize(auxiliaryField, *_normalizer, dim);

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    // :ATTENTION: In quasi-static problems, the time scale is usually quite large
    // (order of tens to hundreds of years), which means that the density scale is very large,
    // and the acceleration scale is very small. Nevertheless, density times gravitational
    // acceleration will have a scale of pressure divided by length and should be within a few orders
    // of magnitude of 1.

    _auxiliaryFactory->addDensity(); // 0
    _auxiliaryFactory->addShearModulus(); // 1
    _auxiliaryFactory->addBulkModulus(); // 2
    if (_gravityField) {
        _auxiliaryFactory->addGravityField(_gravityField);
    } // if
    if (_useBodyForce) {
        _auxiliaryFactory->addBodyForce();
    } // if
    if (_useReferenceState) {
        _auxiliaryFactory->addReferenceStress(); // numA-2
        _auxiliaryFactory->addReferenceStrain(); // numA-1
    } // if

    auxiliaryField->subfieldsSetup();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->zeroLocal();

    assert(_auxiliaryFactory);
    _auxiliaryFactory->initializeSubfields();

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::materials::IsotropicLinearElasticityPlaneStrain::createDerivedField(const pylith::topology::Field& solution,
                                                                            const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;

    PYLITH_METHOD_RETURN(NULL);
} // createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    // Verify solution contains expected fields.
    if (!solution.hasSubfield("displacement")) {
        throw std::runtime_error("Cannot find 'displacement' field in solution; required for material 'IsotropicLinearElasticityPlaneStrain'.");
    } // if
    if (_useInertia && !solution.hasSubfield("velocity")) {
        throw std::runtime_error("Cannot find 'velocity' field in solution; required for material 'IsotropicLinearElasticityPlaneStrain' with inertia.");
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::materials::IsotropicLinearElasticityPlaneStrain::_getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // _getAuxiliaryFactory


// ----------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setKernelsRHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                                                const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsRHSResidual(solution="<<solution.label()<<")");

    std::vector<ResidualKernels> kernels;

    if (!solution.hasSubfield("velocity")) {
        // Displacement
        const PetscPointFunc g0u = (_gravityField && _useBodyForce) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g0v_gravbodyforce :
                                   (_gravityField) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g0v_grav :
                                   (_useBodyForce) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g0v_bodyforce :
                                   NULL;
        const PetscPointFunc g1u = (!_useReferenceState) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g1v : pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g1v_refstate;

        kernels.resize(1);
        kernels[0] = ResidualKernels("displacement", g0u, g1u);
    } else {
        // Displacement
        const PetscPointFunc g0u = pylith::fekernels::DispVel::g0u;
        const PetscPointFunc g1u = NULL;

        // Velocity
        const PetscPointFunc g0v = (_gravityField && _useBodyForce) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g0v_gravbodyforce :
                                   (_gravityField) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g0v_grav :
                                   (_useBodyForce) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g0v_bodyforce :
                                   NULL;
        const PetscPointFunc g1v = (!_useReferenceState) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g1v : pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g1v_refstate;

        kernels.resize(2);
        kernels[0] = ResidualKernels("displacement", g0u, g1u);
        kernels[1] = ResidualKernels("velocity", g0v, g1v);
    } // if/else

    integrator->setKernelsRHSResidual(kernels);

    PYLITH_METHOD_END;
} // _setKernelsRHSResidual


// ----------------------------------------------------------------------
// Set kernels for RHS Jacobian G(t,s).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setKernelsRHSJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                                                const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsRHSJacobian(solution="<<solution.label()<<")");

    std::vector<JacobianKernels> kernels;

    if (!solution.hasSubfield("velocity")) {
        const PetscPointJac Jg0uu = NULL;
        const PetscPointJac Jg1uu = NULL;
        const PetscPointJac Jg2uu = NULL;
        const PetscPointJac Jg3uu = pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jg3vu;

        kernels.resize(1);
        kernels[0] = JacobianKernels("displacement", "displacement", Jg0uu, Jg1uu, Jg2uu, Jg3uu);
    } else {
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
        const PetscPointJac Jg3vu = pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jg3vu;

        const PetscPointJac Jg0vv = NULL;
        const PetscPointJac Jg1vv = NULL;
        const PetscPointJac Jg2vv = NULL;
        const PetscPointJac Jg3vv = NULL;

        kernels.resize(4);
        kernels[0] = JacobianKernels("displacement", "displacement", Jg0uu, Jg1uu, Jg2uu, Jg3uu);
        kernels[1] = JacobianKernels("displacement", "velocity", Jg0uv, Jg1uv, Jg2uv, Jg3uv);
        kernels[2] = JacobianKernels("velocity", "displacement", Jg0vu, Jg1vu, Jg2vu, Jg3vu);
        kernels[3] = JacobianKernels("velocity", "velocity", Jg0vv, Jg1vv, Jg2vv, Jg3vv);

    } // if/else

    integrator->setKernelsRHSJacobian(kernels);

    PYLITH_METHOD_END;
} // setKernelsRHSJacobian


// ----------------------------------------------------------------------
// Set kernels for LHS residual F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setKernelsLHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                                                const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsLHSResidual(solution="<<solution.label()<<")");

    std::vector<ResidualKernels> kernels;

    if (!solution.hasSubfield("velocity")) {
        // F(t,s,\dot{s}) = \vec{0}.
    } else {
        // Displacement
        const PetscPointFunc f0u = pylith::fekernels::DispVel::f0u;
        const PetscPointFunc f1u = NULL;

        // Velocity
        const PetscPointFunc f0v = (_useInertia) ? pylith::fekernels::DispVel::f0v : NULL;
        const PetscPointFunc f1v = NULL;

        kernels.resize(2);
        kernels[0] = ResidualKernels("displacement", f0u, f1u);
        kernels[1] = ResidualKernels("velocity", f0v, f1v);
    } // if/else

    integrator->setKernelsLHSResidual(kernels);

    PYLITH_METHOD_END;
} // setKernelsLHSResidual


// ----------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearElasticityPlaneStrain::_setKernelsLHSJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                                                const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsLHSJacobian(solution="<<solution.label()<<")");

    std::vector<JacobianKernels> kernels;

    if (!solution.hasSubfield("velocity")) {
        const PetscPointJac Jf0uu = pylith::fekernels::DispVel::Jf0uu_zero;
        const PetscPointJac Jf1uu = NULL;
        const PetscPointJac Jf2uu = NULL;
        const PetscPointJac Jf3uu = NULL;

        kernels.resize(1);
        kernels[0] = JacobianKernels("displacement", "displacement", Jf0uu, Jf1uu, Jf2uu, Jf3uu);
    } else {
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

        kernels.resize(4);
        kernels[0] = JacobianKernels("displacement", "displacement", Jf0uu, Jf1uu, Jf2uu, Jf3uu);
        kernels[1] = JacobianKernels("displacement", "velocity", Jf0uv, Jf1uv, Jf2uv, Jf3uv);
        kernels[2] = JacobianKernels("velocity", "displacement", Jf0vu, Jf1vu, Jf2vu, Jf3vu);
        kernels[3] = JacobianKernels("velocity", "velocity", Jf0vv, Jf1vv, Jf2vv, Jf3vv);
    } // if/else

    PYLITH_METHOD_END;
} // setKernelsLHSJacobian


// End of file
