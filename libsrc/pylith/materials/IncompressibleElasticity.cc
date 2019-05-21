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

#include "pylith/materials/IncompressibleElasticity.hh" // implementation of object methods

#include "pylith/materials/RheologyIncompressibleElasticity.hh" // HASA RheologyIncompressibleElasticity
#include "pylith/materials/AuxiliaryFactoryElasticity.hh" // USES AuxiliaryFactoryElasticity
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/fekernels/IncompressibleElasticity.hh" // USES IncompressibleElasticity kernels
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/DispVel.hh" // USES DispVel kernels

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IncompressibleElasticity::IncompressibleElasticity(void) :
    _useBodyForce(false),
    _rheology(NULL) {
    pylith::utils::PyreComponent::setName("incompressibleelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IncompressibleElasticity::~IncompressibleElasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IncompressibleElasticity::deallocate(void) {
    Material::deallocate();

    _rheology = NULL; // :TODO: Use shared pointer.
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
void
pylith::materials::IncompressibleElasticity::useBodyForce(const bool value) {
    PYLITH_COMPONENT_DEBUG("useBodyForce(value="<<value<<")");

    _useBodyForce = value;
} // useBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
bool
pylith::materials::IncompressibleElasticity::useBodyForce(void) const {
    return _useBodyForce;
} // useBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Set bulk rheology.
void
pylith::materials::IncompressibleElasticity::setBulkRheology(pylith::materials::RheologyIncompressibleElasticity* const rheology) {
    PYLITH_COMPONENT_DEBUG("setBulkRheology(rheology="<<rheology<<")");

    _rheology = rheology;
} // setBulkRheology


// ---------------------------------------------------------------------------------------------------------------------
// Get bulk rheology.
pylith::materials::RheologyIncompressibleElasticity*
pylith::materials::IncompressibleElasticity::getBulkRheology(void) const {
    return _rheology;
} // getBulkRheology


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::materials::IncompressibleElasticity::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    // Verify solution contains expected fields.
    if (!solution.hasSubfield("displacement")) {
        throw std::runtime_error("Cannot find 'displacement' field in solution; required for material 'IncompressibleElasticity'.");
    } // if
    if (!solution.hasSubfield("pressure")) {
        throw std::runtime_error("Cannot find 'pressure' field in solution; required for material 'IncompressibleElasticity'.");
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::materials::IncompressibleElasticity::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.label()<<")");

    pylith::feassemble::IntegratorDomain* integrator = new pylith::feassemble::IntegratorDomain(this);assert(integrator);
    integrator->setMaterialId(getMaterialId());

    _setKernelsRHSResidual(integrator, solution);
    _setKernelsRHSJacobian(integrator, solution);
    _setKernelsLHSJacobian(integrator, solution);
    // No state variables.
    // _setKernelsDerivedFields(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::materials::IncompressibleElasticity::createAuxiliaryField(const pylith::topology::Field& solution,
                                                                  const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.label()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->label("IncompressibleElasticity auxiliary field");

    assert(_rheology);
    pylith::materials::AuxiliaryFactoryElasticity* auxiliaryFactory = _rheology->getAuxiliaryFactory();assert(auxiliaryFactory);

    assert(_normalizer);
    auxiliaryFactory->initialize(auxiliaryField, *_normalizer, domainMesh.dimension());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    // :ATTENTION: In quasi-static problems, the time scale is usually quite large
    // (order of tens to hundreds of years), which means that the density scale is very large,
    // and the acceleration scale is very small. Nevertheless, density times gravitational
    // acceleration will have a scale of pressure divided by length and should be within a few orders
    // of magnitude of 1.

    auxiliaryFactory->addDensity(); // 0
    if (_useBodyForce) {
        auxiliaryFactory->addBodyForce();
    } // if
    if (_gravityField) {
        auxiliaryFactory->addGravityField(_gravityField);
    } // if
    _rheology->addAuxiliarySubfields();

    auxiliaryField->subfieldsSetup();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->zeroLocal();

    assert(auxiliaryFactory);
    auxiliaryFactory->setValuesFromDB();

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::materials::IncompressibleElasticity::createDerivedField(const pylith::topology::Field& solution,
                                                                const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.label()<<", domainMesh=)"<<typeid(domainMesh).name()<<") empty method");

    PYLITH_METHOD_RETURN(NULL);
} // createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::materials::IncompressibleElasticity::_getAuxiliaryFactory(void) {
    return _rheology->getAuxiliaryFactory();
} // _getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::materials::IncompressibleElasticity::_setKernelsRHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                                    const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsRHSResidual(integrator="<<integrator<<", solution="<<solution.label()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.mesh().coordsys();

    std::vector<ResidualKernels> kernels(2);

    // Displacement
    const PetscPointFunc g0u =
        (_gravityField && _useBodyForce) ? pylith::fekernels::Elasticity::g0v_gravbodyforce :
        (_gravityField) ? pylith::fekernels::Elasticity::g0v_grav :
        (_useBodyForce) ? pylith::fekernels::Elasticity::g0v_bodyforce :
        NULL;
    const PetscPointFunc g1u = _rheology->getKernelRHSResidualStress(coordsys);

    // Pressure
    const PetscPointFunc g0p = _rheology->getKernelRHSResidualPressure(coordsys);
    const PetscPointFunc g1p = NULL;

    kernels[0] = ResidualKernels("displacement", g0u, g1u);
    kernels[1] = ResidualKernels("pressure", g0p, g1p);

    assert(integrator);
    integrator->setKernelsRHSResidual(kernels);

    PYLITH_METHOD_END;
} // _setKernelsRHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS Jacobian G(t,s).
void
pylith::materials::IncompressibleElasticity::_setKernelsRHSJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                                    const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsRHSJacobian(integrator="<<integrator<<", solution="<<solution.label()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.mesh().coordsys();

    std::vector<JacobianKernels> kernels(4);

    const PetscPointJac Jg0uu = NULL;
    const PetscPointJac Jg1uu = NULL;
    const PetscPointJac Jg2uu = NULL;
    const PetscPointJac Jg3uu = _rheology->getKernelRHSJacobianElasticConstants(coordsys);

    const PetscPointJac Jg0up = NULL;
    const PetscPointJac Jg1up = NULL;
    const PetscPointJac Jg2up = pylith::fekernels::IncompressibleElasticity::Jg2up;
    const PetscPointJac Jg3up = NULL;

    const PetscPointJac Jg0pu = NULL;
    const PetscPointJac Jg1pu = pylith::fekernels::IncompressibleElasticity::Jg1pu;
    const PetscPointJac Jg2pu = NULL;
    const PetscPointJac Jg3pu = NULL;

    const PetscPointJac Jg0pp = _rheology->getKernelRHSJacobianInverseBulkModulus(coordsys);
    const PetscPointJac Jg1pp = NULL;
    const PetscPointJac Jg2pp = NULL;
    const PetscPointJac Jg3pp = NULL;

    kernels[0] = JacobianKernels("displacement", "displacement", Jg0uu, Jg1uu, Jg2uu, Jg3uu);
    kernels[1] = JacobianKernels("displacement", "pressure", Jg0up, Jg1up, Jg2up, Jg3up);
    kernels[2] = JacobianKernels("pressure", "displacement", Jg0pu, Jg1pu, Jg2pu, Jg3pu);
    kernels[3] = JacobianKernels("pressure", "pressure", Jg0pp, Jg1pp, Jg2pp, Jg3pp);

    assert(integrator);
    integrator->setKernelsRHSJacobian(kernels);

    PYLITH_METHOD_END;
} // setKernelsRHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,s,\dot{s}).
void
pylith::materials::IncompressibleElasticity::_setKernelsLHSJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                                    const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsLHSJacobian(integrator="<<integrator<<", solution="<<solution.label()<<")");

    std::vector<JacobianKernels> kernels;

    const PetscPointJac Jf0 = pylith::fekernels::DispVel::Jf0uu_zero;
    const PetscPointJac Jf1 = NULL;
    const PetscPointJac Jf2 = NULL;
    const PetscPointJac Jf3 = NULL;

    kernels.resize(4);
    kernels[0] = JacobianKernels("displacement", "displacement", Jf0, Jf1, Jf2, Jf3);
    kernels[1] = JacobianKernels("displacement", "pressure", Jf0, Jf1, Jf2, Jf3);
    kernels[2] = JacobianKernels("pressure", "displacement", Jf0, Jf1, Jf2, Jf3);
    kernels[3] = JacobianKernels("pressure", "pressure", Jf0, Jf1, Jf2, Jf3);

    assert(integrator);
    integrator->setKernelsLHSJacobian(kernels);

    PYLITH_METHOD_END;
} // setKernelsLHSJacobian


// End of file
