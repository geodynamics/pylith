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

#include "pylith/materials/IsotropicLinearIncompElasticity.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryElastic.hh" // USES AuxiliaryFactory
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
//#include "pylith/fekernels/IsotropicLinearIncompElasticity.hh" // USES IsotropicLinearIncompElasticity kernels
#include "pylith/fekernels/DispVel.hh" // USES DispVel kernels

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearIncompElasticity::IsotropicLinearIncompElasticity(void) :
    _useInertia(false),
    _useBodyForce(false),
    _useReferenceState(false) {
    pylith::utils::PyreComponent::name("isotopiclinearincomplinearelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearIncompElasticity::~IsotropicLinearIncompElasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearIncompElasticity::deallocate(void) {
    Material::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Include inertia?
void
pylith::materials::IsotropicLinearIncompElasticity::useInertia(const bool value) {
    PYLITH_COMPONENT_DEBUG("useInertia(value="<<value<<")");
    _useInertia = value;
} // useInertia


// ---------------------------------------------------------------------------------------------------------------------
// Include inertia?
bool
pylith::materials::IsotropicLinearIncompElasticity::useInertia(void) const {
    return _useInertia;
} // useInertia


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
void
pylith::materials::IsotropicLinearIncompElasticity::useBodyForce(const bool value) {
    PYLITH_COMPONENT_DEBUG("useBodyForce(value="<<value<<")");

    _useBodyForce = value;
} // useBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
bool
pylith::materials::IsotropicLinearIncompElasticity::useBodyForce(void) const {
    return _useBodyForce;
} // useBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
void
pylith::materials::IsotropicLinearIncompElasticity::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
bool
pylith::materials::IsotropicLinearIncompElasticity::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::materials::IsotropicLinearIncompElasticity::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    // Verify solution contains expected fields.
    if (!solution.hasSubfield("displacement")) {
        throw std::runtime_error("Cannot find 'displacement' field in solution; required for material 'IsotropicLinearIncompElasticity'.");
    } // if
    if (!solution.hasSubfield("pressure")) {
        throw std::runtime_error("Cannot find 'pressure' field in solution; required for material 'IsotropicLinearIncompElasticity'.");
    } // if
    if (_useInertia && !solution.hasSubfield("velocity")) {
        throw std::runtime_error("Cannot find 'velocity' field in solution; required for material 'IsotropicLinearIncompElasticity' with inertia.");
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::materials::IsotropicLinearIncompElasticity::createIntegrator(const pylith::topology::Field& solution) {
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
pylith::materials::IsotropicLinearIncompElasticity::createAuxiliaryField(const pylith::topology::Field& solution,
                                                                         const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->label("IsotropicLinearIncompElasticity auxiliary field");

    assert(_auxiliaryFactory);
    assert(_normalizer);
    _auxiliaryFactory->initialize(auxiliaryField, *_normalizer, domainMesh.dimension());

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
pylith::materials::IsotropicLinearIncompElasticity::createDerivedField(const pylith::topology::Field& solution,
                                                                       const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;

    PYLITH_METHOD_RETURN(NULL);
} // createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::materials::IsotropicLinearIncompElasticity::_getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // _getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::materials::IsotropicLinearIncompElasticity::_setKernelsRHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                                           const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsRHSResidual(integrator="<<integrator<<", solution="<<solution.label()<<")");

    //const int spaceDim = solution.spaceDim();

    std::vector<ResidualKernels> kernels;

#if 1
    PYLITH_COMPONENT_ERROR(":TODO: @charles Implement gravbodyforce, grav, and bodyforce kernels.");
#else
    if (!solution.hasSubfield("velocity")) {
        // Displacement
        const PetscPointFunc g0u = (_gravityField && _useBodyForce) ? pylith::fekernels::IsotropicLinearIncompElasticity::g0v_gravbodyforce :
                                   (_gravityField) ? pylith::fekernels::IsotropicLinearIncompElasticity::g0v_grav :
                                   (_useBodyForce) ? pylith::fekernels::IsotropicLinearIncompElasticity::g0v_bodyforce :
                                   NULL;
        const PetscPointFunc g1u = (!_useReferenceState) ? pylith::fekernels::IsotropicLinearIncompElasticity::g1v : pylith::fekernels::IsotropicLinearIncompElasticity::g1v_refstate;

        // Pressure
        const PetscPointFunc g0p = pylith::fekernels::IncompressibleElasticity::g0p;
        const PetscPointFunc g1p = NULL;

        kernels.resize(1);
        kernels[0] = ResidualKernels("displacement", g0u, g1u);
        kernels[1] = ResidualKernels("pressure", g0p, g1p);

    } else {
        const PetscInt i_vel = solution.subfieldInfo("velocity").index;

        // Displacement
        const PetscPointFunc g0u = pylith::fekernels::DispVel::g0u;
        const PetscPointFunc g1u = NULL;

        // Pressure
        const PetscPointFunc g0p = pylith::fekernels::IncompressibleElasticity::g0p;
        const PetscPointFunc g1p = NULL;

        // Velocity
        const PetscPointFunc g0v = (_gravityField && _useBodyForce) ? pylith::fekernels::IsotropicLinearIncompElasticity::g0v_gravbodyforce :
                                   (_gravityField) ? pylith::fekernels::IsotropicLinearIncompElasticity::g0v_grav :
                                   (_useBodyForce) ? pylith::fekernels::IsotropicLinearIncompElasticity::g0v_bodyforce :
                                   NULL;
        const PetscPointFunc g1v = (!_useReferenceState) ? pylith::fekernels::IsotropicLinearIncompElasticity::g1v : pylith::fekernels::IsotropicLinearIncompElasticity::g1v_refstate;

        kernels.resize(3);
        kernels[0] = ResidualKernels("displacement", g0u, g1u);
        kernels[1] = ResidualKernels("velocity", g0v, g1v);
        kernels[2] = ResidualKernels("pressure", g0p, g1p);
    } // if/else
#endif

    assert(integrator);
    integrator->setKernelsRHSResidual(kernels);

    PYLITH_METHOD_END;
} // _setKernelsRHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS Jacobian G(t,s).
void
pylith::materials::IsotropicLinearIncompElasticity::_setKernelsRHSJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                                           const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsRHSJacobian(integrator="<<integrator<<", solution="<<solution.label()<<")");

    //const int spaceDim = solution.spaceDim();

    std::vector<JacobianKernels> kernels;

#if 1
    PYLITH_COMPONENT_ERROR(":TODO: @charles Implement RHS Jacobian kernels.");
#else
    if (!solution.hasSubfield("velocity")) {
        const PetscPointJac Jg0uu = NULL;
        const PetscPointJac Jg1uu = NULL;
        const PetscPointJac Jg2uu = NULL;
        const PetscPointJac Jg3uu = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::Jg3vu :
                                    (2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jg3vu :
                                    NULL;

        const PetscPointJac Jg0up = NULL;
        const PetscPointJac Jg1up = NULL;
        const PetscPointJac Jg2up = NULL;
        const PetscPointJac Jg3up = NULL;

        const PetscPointJac Jg0pu = NULL;
        const PetscPointJac Jg1pu = pylith::fekernels::IsotropicLinearIncompElasticity::Jg1pu;
        const PetscPointJac Jg2pu = NULL;
        const PetscPointJac Jg3pu = NULL;

        const PetscPointJac Jg0pp = pylith::fekernels::IncompressibleElasticity::Jg0pp;
        const PetscPointJac Jg1pp = NULL;
        const PetscPointJac Jg2pp = NULL;
        const PetscPointJac Jg3pp = NULL;

        kernels.resize(4);
        kernels[0] = JacobianKernels("displacement", "displacement", Jg0uu, Jg1uu, Jg2uu, Jg3uu);
        kernels[1] = JacobianKernels("displacement", "pressure", Jg0up, Jg1up, Jg2up, Jg3up);
        kernels[2] = JacobianKernels("pressure", "displacement", Jg0pu, Jg1pu, Jg2pu, Jg3pu);
        kernels[3] = JacobianKernels("pressure", "pressure", Jg0pp, Jg1pp, Jg2pp, Jg3pp);

    } else {
        PYLITH_COMPONENT_ERROR("IsotropicLinearIncompElasticity with velocity solution field not implemented.");
        throw std::logic_error("IsotropicLinearIncompElasticity with velocity solution field not implemented.");
    } // if/else
#endif

    PYLITH_METHOD_END;
} // _setKernelsRHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for LHS residual F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearIncompElasticity::_setKernelsLHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                                           const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsLHSResidual(integrator="<<integrator<<", solution="<<solution.label()<<")");

    std::vector<ResidualKernels> kernels;

    if (!solution.hasSubfield("velocity")) {
        // F(t,s,\dot{s}) = \vec{0}.
    } else {
        PYLITH_COMPONENT_ERROR("IsotropicLinearIncompElasticity with velocity solution field not implemented.");
        throw std::logic_error("IsotropicLinearIncompElasticity with velocity solution field not implemented.");
    } // if/else

    assert(integrator);
    integrator->setKernelsLHSResidual(kernels);

    PYLITH_METHOD_END;
} // _setKernelsLHSResidual


// ----------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearIncompElasticity::_setKernelsLHSJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                                           const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsLHSJacobianImplicit(integrator="<<integrator<<", solution="<<solution.label()<<")");

    std::vector<JacobianKernels> kernels;

    if (!solution.hasSubfield("velocity")) {
        const PetscPointJac Jf0uu = pylith::fekernels::DispVel::Jf0uu_zero;
        const PetscPointJac Jf1uu = NULL;
        const PetscPointJac Jf2uu = NULL;
        const PetscPointJac Jf3uu = NULL;

        const PetscPointJac Jf0up = NULL;
        const PetscPointJac Jf1up = NULL;
        const PetscPointJac Jf2up = NULL;
        const PetscPointJac Jf3up = NULL;

        const PetscPointJac Jf0pu = NULL;
        const PetscPointJac Jf1pu = NULL;
        const PetscPointJac Jf2pu = NULL;
        const PetscPointJac Jf3pu = NULL;

        const PetscPointJac Jf0pp = NULL;
        const PetscPointJac Jf1pp = NULL;
        const PetscPointJac Jf2pp = NULL;
        const PetscPointJac Jf3pp = NULL;

        kernels.resize(4);
        kernels[0] = JacobianKernels("displacement", "displacement", Jf0uu, Jf1uu, Jf2uu, Jf3uu);
        kernels[1] = JacobianKernels("displacement", "pressure", Jf0up, Jf1up, Jf2up, Jf3up);
        kernels[2] = JacobianKernels("pressure", "displacement", Jf0pu, Jf1pu, Jf2pu, Jf3pu);
        kernels[3] = JacobianKernels("pressure", "pressure", Jf0pp, Jf1pp, Jf2pp, Jf3pp);
    } else {
        PYLITH_COMPONENT_ERROR("IsotropicLinearIncompElasticity with velocity solution field not implemented.");
        throw std::logic_error("IsotropicLinearIncompElasticity with velocity solution field not implemented.");
    } // if/else

    PYLITH_METHOD_END;
} // _setKernelsLHSJacobian


// End of file
