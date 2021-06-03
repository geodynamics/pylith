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

#include "pylith/materials/Poroelasticity.hh" // implementation of object methods

#include "pylith/materials/RheologyPoroelasticity.hh" // HASA RheologyPoroelasticity
#include "pylith/materials/AuxiliaryFactoryPoroelastic.hh" // USES AuxiliaryFactory
#include "pylith/materials/DerivedFactoryElasticity.hh" // USES DerivedFactoryElasticity
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/fekernels/Poroelasticity.hh" // USES Poroelasticity kernels
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/DispVel.hh" // USES DispVel kernels

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::Poroelasticity::Poroelasticity(void) :
    _useInertia(false),
    _useBodyForce(false),
    _useReferenceState(false),
    _useSourceDensity(false),
    _rheology(NULL),
    _derivedFactory(new pylith::materials::DerivedFactoryElasticity) {
    pylith::utils::PyreComponent::setName("poroelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::Poroelasticity::~Poroelasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::Poroelasticity::deallocate(void) {
    Material::deallocate();

    delete _derivedFactory;_derivedFactory = NULL;
    _rheology = NULL; // :TODO: Use shared pointer.
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
void
pylith::materials::Poroelasticity::useBodyForce(const bool value) {
    PYLITH_COMPONENT_DEBUG("useBodyForce(value="<<value<<")");

    _useBodyForce = value;
} // useBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
bool
pylith::materials::Poroelasticity::useBodyForce(void) const {
    return _useBodyForce;
} // useBodyForce


// ----------------------------------------------------------------------
// Include source density?
void
pylith::materials::Poroelasticity::useSourceDensity(const bool value) {
    PYLITH_COMPONENT_DEBUG("useSourceDensity(value="<<value<<")");

    _useSourceDensity = value;
} // useSourceDensity


// ----------------------------------------------------------------------
// Include source density?
bool
pylith::materials::Poroelasticity::useSourceDensity(void) const {
    return _useSourceDensity;
} // useSourceDensity


// ---------------------------------------------------------------------------------------------------------------------
// Set bulk rheology.
void
pylith::materials::Poroelasticity::setBulkRheology(pylith::materials::RheologyPoroelasticity* const rheology) {
    PYLITH_COMPONENT_DEBUG("setBulkRheology(rheology="<<rheology<<")");

    _rheology = rheology;
} // setBulkRheology


// ---------------------------------------------------------------------------------------------------------------------
// Get bulk rheology.
pylith::materials::RheologyPoroelasticity*
pylith::materials::Poroelasticity::getBulkRheology(void) const {
    return _rheology;
} // getBulkRheology


// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::materials::Poroelasticity::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    // Verify solution contains expected fields.
    if (!solution.hasSubfield("displacement")) {
        throw std::runtime_error("Cannot find 'displacement' field in solution; required for material 'Poroelasticity'.");
    } // if
    if (!solution.hasSubfield("pressure")) {
        throw std::runtime_error("Cannot find 'pressure' field in solution; required for material 'Poroelasticity'.");
    } // if
    switch (_formulation) {
    case QUASISTATIC:
        if (!_useInertia && !solution.hasSubfield("trace_strain")) {
            throw std::runtime_error("Cannot find 'trace_strain' field in solution; required for material 'Poroelasticity'.");
        } // if
        break;
    case DYNAMIC:
    case DYNAMIC_IMEX:
        if (_useInertia && !solution.hasSubfield("velocity")) {
            throw std::runtime_error("Cannot find 'velocity' field in solution; required for material 'Poroelasticity' with inertia.");
        } // if
    } // switch
    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::materials::Poroelasticity::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    pylith::feassemble::IntegratorDomain* integrator = new pylith::feassemble::IntegratorDomain(this);assert(integrator);
    integrator->setLabelName("material-id");
    integrator->setLabelValue(getMaterialId());

    switch (_formulation) {
    case QUASISTATIC:
        _useInertia = false;
        break;
    case DYNAMIC:
        _useInertia = true;
        break;
    case DYNAMIC_IMEX:
        _useInertia = true;
        integrator->setLHSJacobianTriggers(pylith::feassemble::Integrator::NEW_JACOBIAN_TIME_STEP_CHANGE);
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    _setKernelsLHSResidual(integrator, solution);
    _setKernelsRHSResidual(integrator, solution);
    _setKernelsLHSJacobian(integrator, solution);
    _setKernelsUpdateStateVars(integrator, solution);
    _setKernelsDerivedField(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::materials::Poroelasticity::createAuxiliaryField(const pylith::topology::Field& solution,
                                                        const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->setLabel("Poroelasticity auxiliary field");

    assert(_rheology);
    pylith::materials::AuxiliaryFactoryPoroelasticity* auxiliaryFactory = _rheology->getAuxiliaryFactory();assert(auxiliaryFactory);

    assert(_normalizer);
    auxiliaryFactory->initialize(auxiliaryField, *_normalizer, domainMesh.dimension());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    // :ATTENTION: In quasi-static problems, the time scale is usually quite large
    // (order of tens to hundreds of years), which means that the density scale is very large,
    // and the acceleration scale is very small. Nevertheless, density times gravitational
    // acceleration will have a scale of pressure divided by length and should be within a few orders
    // of magnitude of 1.

    // ---------------------------------
    // Required Auxiliary
    auxiliaryFactory->addSolidDensity(); // 0 Rock Density
    auxiliaryFactory->addFluidDensity(); // 1 Fluid Density
    auxiliaryFactory->addFluidViscosity(); // 2 Fluid Viscosity
    auxiliaryFactory->addPorosity(); // 3 Porosity

    // ---------------------------------
    // Optional Auxiliary
    if (_useBodyForce) {
        auxiliaryFactory->addBodyForce(); // +1
    } // if
    if (_gravityField) {
        auxiliaryFactory->addGravityField(_gravityField); // +1
    } // if
    if (_useSourceDensity) {
        auxiliaryFactory->addSourceDensity(); // +1
    } // if
    _rheology->addAuxiliarySubfields();

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->zeroLocal();
    auxiliaryField->createOutputVector();

    assert(auxiliaryFactory);
    auxiliaryFactory->setValuesFromDB();

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::materials::Poroelasticity::createDerivedField(const pylith::topology::Field& solution,
                                                      const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    assert(_derivedFactory);
    if (_derivedFactory->getNumSubfields() == 1) {
        PYLITH_METHOD_RETURN(NULL);
    } // if

    pylith::topology::Field* derivedField = new pylith::topology::Field(domainMesh);assert(derivedField);
    derivedField->setLabel("Poroelasticity derived field");

    assert(_normalizer);
    _derivedFactory->initialize(derivedField, *_normalizer, domainMesh.dimension());
    _derivedFactory->addSubfields();

    derivedField->subfieldsSetup();
    derivedField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *derivedField);
    derivedField->allocate();
    derivedField->zeroLocal();
    derivedField->createOutputVector();

    PYLITH_METHOD_RETURN(derivedField);
} // createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::materials::Poroelasticity::_getAuxiliaryFactory(void) {
    assert(_rheology);
    return _rheology->getAuxiliaryFactory();
} // _getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::Poroelasticity::_updateKernelConstants(const PylithReal dt) {
    assert(_rheology);
    _rheology->updateKernelConstants(&_kernelConstants, dt);
} // _updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Get derived factory associated with physics.
pylith::topology::FieldFactory*
pylith::materials::Poroelasticity::_getDerivedFactory(void) {
    return _derivedFactory;
} // _getDerivedFactory


// ----------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
// RHS residual is for explicit time stepping (dynamic).
void
pylith::materials::Poroelasticity::_setKernelsRHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                          const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsRHSResidual(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.mesh().getCoordSys();

    std::vector<ResidualKernels> kernels;

    switch (_formulation) {
    case QUASISTATIC: {
        break;
    } // QUASISTATIC

    case DYNAMIC_IMEX:
    case DYNAMIC: {
        // Displacement
        const PetscPointFunc g0u = pylith::fekernels::Poroelasticity::g0u;
        const PetscPointFunc g1u = NULL;

        // Pressure
        const PetscPointFunc g0p = _rheology->getKernelg0p(coordsys, _useBodyForce, _gravityField, _useSourceDensity);
        const PetscPointFunc g1p = _rheology->getKernelg1p_explicit(coordsys, _gravityField); // darcy velocity

        // Velocity
        PetscPointFunc g0v = NULL;
        const PetscPointFunc g1v = _rheology->getKernelg1v_explicit(coordsys);

        const int bitBodyForce = _useBodyForce ? 0x1 : 0x0;
        const int bitGravity = _gravityField ? 0x2 : 0x0;
        const int bitSourceDensity = _useSourceDensity ? 0x4 : 0x0;
        const int bitUse = bitBodyForce | bitGravity | bitSourceDensity;

        switch (bitUse) {
        case 0x0:
            break;
        case 0x1:
            g0v = pylith::fekernels::Poroelasticity::g0v_bodyforce;
            break;
        case 0x2:
            g0v = pylith::fekernels::Poroelasticity::g0v_grav;
            break;
        case 0x4:
            break;
        case 0x3:
            g0v = pylith::fekernels::Poroelasticity::g0v_grav_bodyforce;
            break;
        case 0x5:
            g0v = pylith::fekernels::Poroelasticity::g0v_bodyforce;
            break;
        case 0x6:
            g0v = pylith::fekernels::Poroelasticity::g0v_grav;
            break;
        case 0x7:
            g0v = pylith::fekernels::Poroelasticity::g0v_grav_bodyforce;
            break;
        default:
            PYLITH_COMPONENT_LOGICERROR("Unknown case (bitUse=" << bitUse << ") for Poroelasticity RHS residual kernels.");
        } // switch

        kernels.resize(3);
        kernels[0] = ResidualKernels("displacement", g0u, g1u);
        kernels[1] = ResidualKernels("pressure", g0p, g1p);
        kernels[2] = ResidualKernels("velocity", g0v, g1v);
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    assert(integrator);
    integrator->setKernelsRHSResidual(kernels);

    PYLITH_METHOD_END;
} // _setKernelsRHSResidual


// ----------------------------------------------------------------------
// Set kernels for LHS residual F(t,s,\dot{s}).
void
pylith::materials::Poroelasticity::_setKernelsLHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                          const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsLHSResidual(integrator="<<integrator<<",solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.mesh().getCoordSys();

    std::vector<ResidualKernels> kernels;
    switch (_formulation) {
    case QUASISTATIC: {
        // Displacement
        PetscPointFunc f0u = NULL;
        const int bitBodyForce = _useBodyForce ? 0x1 : 0x0;
        const int bitGravity = _gravityField ? 0x2 : 0x0;
        const int bitSourceDensity = _useSourceDensity ? 0x4 : 0x0;
        const int bitUse = bitBodyForce | bitGravity | bitSourceDensity;

        switch (bitUse) {
        case 0x0:
            break;
        case 0x1:
            f0u = pylith::fekernels::Poroelasticity::g0v_bodyforce;
            break;
        case 0x2:
            f0u = pylith::fekernels::Poroelasticity::g0v_grav;
            break;
        case 0x4:
            break;
        case 0x3:
            f0u = pylith::fekernels::Poroelasticity::g0v_grav_bodyforce;
            break;
        case 0x5:
            f0u = pylith::fekernels::Poroelasticity::g0v_bodyforce;
            break;
        case 0x6:
            f0u = pylith::fekernels::Poroelasticity::g0v_grav;
            break;
        case 0x7:
            f0u = pylith::fekernels::Poroelasticity::g0v_grav_bodyforce;
            break;
        default:
            PYLITH_COMPONENT_LOGICERROR("Unknown case (bitUse=" << bitUse << ") for Poroelasticity RHS residual kernels.");
        } // switch
        const PetscPointFunc f1u = _rheology->getKernelf1u_implicit(coordsys);

        // Pressure
        PetscPointFunc f0p = _rheology->getKernelf0p_implicit(coordsys, _useBodyForce, _gravityField, _useSourceDensity);
        PetscPointFunc f1p = _rheology->getKernelf1p_implicit(coordsys, _gravityField); // darcy velocity

        // Volumetric Strain
        const PetscPointFunc f0e = pylith::fekernels::Poroelasticity::f0e;
        const PetscPointFunc f1e = NULL;

        kernels.resize(3);
        kernels[0] = ResidualKernels("displacement", f0u, f1u);
        kernels[1] = ResidualKernels("pressure",     f0p, f1p);
        kernels[2] = ResidualKernels("trace_strain", f0e, f1e);
        break;
    } // QUASISTATIC

    case DYNAMIC_IMEX: {
        // Displacement
        const PetscPointFunc f0u = pylith::fekernels::DispVel::f0u;
        const PetscPointFunc f1u = NULL;

        // Pressure
        const PetscPointFunc f0p = _rheology->getKernelf0p_explicit(coordsys);
        const PetscPointFunc f1p = NULL;

        // Velocity
        assert(_useInertia);
        const PetscPointFunc f0v = pylith::fekernels::Poroelasticity::f0v;
        const PetscPointFunc f1v = NULL;

        kernels.resize(3);
        kernels[0] = ResidualKernels("displacement", f0u, f1u);
        kernels[1] = ResidualKernels("pressure",     f0p, f1p);
        kernels[2] = ResidualKernels("velocity",     f0v, f1v);
        break;
    } // DYNAMIC_IMEX
    case DYNAMIC: {
        // Displacement
        const PetscPointFunc f0u = pylith::fekernels::DispVel::f0u;
        const PetscPointFunc f1u = NULL;

        // Pressure
        const PetscPointFunc f0p = _rheology->getKernelf0p_explicit(coordsys);
        const PetscPointFunc f1p = NULL;

        // Velocity
        assert(_useInertia);
        const PetscPointFunc f0v = pylith::fekernels::Poroelasticity::f0v;
        const PetscPointFunc f1v = NULL;

        kernels.resize(3);
        kernels[0] = ResidualKernels("displacement", f0u, f1u);
        kernels[1] = ResidualKernels("pressure",     f0p, f1p);
        kernels[2] = ResidualKernels("velocity",     f0v, f1v);
        break;
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    assert(integrator);
    integrator->setKernelsLHSResidual(kernels);

    PYLITH_METHOD_END;
} // _setKernelsLHSResidual


// ----------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,s,\dot{s}).
void
pylith::materials::Poroelasticity::_setKernelsLHSJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                          const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsLHSJacobian(integrator="<<integrator<<",solution="<<solution.getLabel()<<")");
    const spatialdata::geocoords::CoordSys* coordsys = solution.mesh().getCoordSys();
    std::vector<JacobianKernels> kernels(9);

    switch (_formulation) {
    case QUASISTATIC: {
        const PetscPointJac Jf0uu = NULL;
        const PetscPointJac Jf1uu = NULL;
        const PetscPointJac Jf2uu = NULL;
        const PetscPointJac Jf3uu = _rheology->getKernelJf3uu(coordsys);

        const PetscPointJac Jf0up = NULL;
        const PetscPointJac Jf1up = NULL;
        const PetscPointJac Jf2up = _rheology->getKernelJf2up(coordsys);
        const PetscPointJac Jf3up = NULL;

        const PetscPointJac Jf0ue = NULL;
        const PetscPointJac Jf1ue = NULL;
        const PetscPointJac Jf2ue = _rheology->getKernelJf2ue(coordsys);
        const PetscPointJac Jf3ue = NULL;

        const PetscPointJac Jf0pu = NULL;
        const PetscPointJac Jf1pu = NULL;
        const PetscPointJac Jf2pu = NULL;
        const PetscPointJac Jf3pu = NULL;

        const PetscPointJac Jf0pp = _rheology->getKernelJf0pp(coordsys);
        const PetscPointJac Jf1pp = NULL;
        const PetscPointJac Jf2pp = NULL;
        const PetscPointJac Jf3pp = _rheology->getKernelJf3pp(coordsys);

        const PetscPointJac Jf0pe = _rheology->getKernelJf0pe(coordsys);
        const PetscPointJac Jf1pe = NULL;
        const PetscPointJac Jf2pe = NULL;
        const PetscPointJac Jf3pe = NULL;

        const PetscPointJac Jf0eu = NULL;
        const PetscPointJac Jf1eu = pylith::fekernels::Poroelasticity::Jf1eu;
        const PetscPointJac Jf2eu = NULL;
        const PetscPointJac Jf3eu = NULL;

        const PetscPointJac Jf0ep = NULL;
        const PetscPointJac Jf1ep = NULL;
        const PetscPointJac Jf2ep = NULL;
        const PetscPointJac Jf3ep = NULL;

        const PetscPointJac Jf0ee = pylith::fekernels::Poroelasticity::Jf0ee;
        const PetscPointJac Jf1ee = NULL;
        const PetscPointJac Jf2ee = NULL;
        const PetscPointJac Jf3ee = NULL;

        kernels[0] = JacobianKernels("displacement",  "displacement",  Jf0uu, Jf1uu, Jf2uu, Jf3uu);
        kernels[1] = JacobianKernels("displacement",  "pressure",      Jf0up, Jf1up, Jf2up, Jf3up);
        kernels[2] = JacobianKernels("displacement",  "trace_strain",  Jf0ue, Jf1ue, Jf2ue, Jf3ue);
        kernels[3] = JacobianKernels("pressure",      "displacement",  Jf0pu, Jf1pu, Jf2pu, Jf3pu);
        kernels[4] = JacobianKernels("pressure",      "pressure",      Jf0pp, Jf1pp, Jf2pp, Jf3pp);
        kernels[5] = JacobianKernels("pressure",      "trace_strain",  Jf0pe, Jf1pe, Jf2pe, Jf3pe);
        kernels[6] = JacobianKernels("trace_strain",  "displacement",  Jf0eu, Jf1eu, Jf2eu, Jf3eu);
        kernels[7] = JacobianKernels("trace_strain",  "pressure",      Jf0ep, Jf1ep, Jf2ep, Jf3ep);
        kernels[8] = JacobianKernels("trace_strain",  "trace_strain",  Jf0ee, Jf1ee, Jf2ee, Jf3ee);
        break;
    } // QUASISTATIC
    case DYNAMIC_IMEX:
    case DYNAMIC: {
        const PetscPointJac Jf0uu = pylith::fekernels::DispVel::Jf0uu_stshift;
        const PetscPointJac Jf1uu = NULL;
        const PetscPointJac Jf2uu = NULL;
        const PetscPointJac Jf3uu = NULL;

        const PetscPointJac Jf0up = NULL;
        const PetscPointJac Jf1up = NULL;
        const PetscPointJac Jf2up = NULL;
        const PetscPointJac Jf3up = NULL;

        const PetscPointJac Jf0uv = NULL;
        const PetscPointJac Jf1uv = NULL;
        const PetscPointJac Jf2uv = NULL;
        const PetscPointJac Jf3uv = NULL;

        const PetscPointJac Jf0pu = NULL;
        const PetscPointJac Jf1pu = NULL;
        const PetscPointJac Jf2pu = NULL;
        const PetscPointJac Jf3pu = NULL;

        const PetscPointJac Jf0pp = _rheology->getKernelJf0pp(coordsys);
        const PetscPointJac Jf1pp = NULL;
        const PetscPointJac Jf2pp = NULL;
        const PetscPointJac Jf3pp = NULL;

        const PetscPointJac Jf0pv = NULL;
        const PetscPointJac Jf1pv = NULL;
        const PetscPointJac Jf2pv = NULL;
        const PetscPointJac Jf3pv = NULL;

        const PetscPointJac Jf0vu = NULL;
        const PetscPointJac Jf1vu = NULL;
        const PetscPointJac Jf2vu = NULL;
        const PetscPointJac Jf3vu = NULL;

        const PetscPointJac Jf0vp = NULL;
        const PetscPointJac Jf1vp = NULL;
        const PetscPointJac Jf2vp = NULL;
        const PetscPointJac Jf3vp = NULL;

        const PetscPointJac Jf0vv = pylith::fekernels::Poroelasticity::Jf0vv;
        const PetscPointJac Jf1vv = NULL;
        const PetscPointJac Jf2vv = NULL;
        const PetscPointJac Jf3vv = NULL;

        kernels[0] = JacobianKernels("displacement",  "displacement",  Jf0uu, Jf1uu, Jf2uu, Jf3uu);
        kernels[1] = JacobianKernels("displacement",  "pressure",      Jf0up, Jf1up, Jf2up, Jf3up);
        kernels[2] = JacobianKernels("displacement",  "velocity",      Jf0uv, Jf1uv, Jf2uv, Jf3uv);
        kernels[3] = JacobianKernels("pressure",      "displacement",  Jf0pu, Jf1pu, Jf2pu, Jf3pu);
        kernels[4] = JacobianKernels("pressure",      "pressure",      Jf0pp, Jf1pp, Jf2pp, Jf3pp);
        kernels[5] = JacobianKernels("pressure",      "velocity",      Jf0pv, Jf1pv, Jf2pv, Jf3pv);
        kernels[6] = JacobianKernels("velocity",      "displacement",  Jf0vu, Jf1vu, Jf2vu, Jf3vu);
        kernels[7] = JacobianKernels("velocity",      "pressure",      Jf0vp, Jf1vp, Jf2vp, Jf3vp);
        kernels[8] = JacobianKernels("velocity",      "velocity",      Jf0vv, Jf1vv, Jf2vv, Jf3vv);
        break;
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    assert(integrator);
    integrator->setKernelsLHSJacobian(kernels);

    PYLITH_METHOD_END;
} // _setKernelsLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for computing updated state variables in auxiliary field.
void
pylith::materials::Poroelasticity::_setKernelsUpdateStateVars(pylith::feassemble::IntegratorDomain* integrator,
                                                              const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsUpdateStateVars(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.mesh().getCoordSys();
    assert(coordsys);

    std::vector<ProjectKernels> kernels;
    _rheology->addKernelsUpdateStateVars(&kernels, coordsys);

    integrator->setKernelsUpdateStateVars(kernels);

    PYLITH_METHOD_END;
} // _setKernelsUpdateStateVars


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for computing derived field.
void
pylith::materials::Poroelasticity::_setKernelsDerivedField(pylith::feassemble::IntegratorDomain* integrator,
                                                           const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsDerivedField(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.mesh().getCoordSys();
    assert(coordsys);

    std::vector<ProjectKernels> kernels(2);
    kernels[0] = ProjectKernels("cauchy_stress", _rheology->getKernelDerivedCauchyStress(coordsys));

    const int spaceDim = coordsys->getSpaceDim();
    const PetscPointFunc strainKernel =
        (3 == spaceDim) ? pylith::fekernels::Poroelasticity3D::cauchyStrain :
        (2 == spaceDim) ? pylith::fekernels::PoroelasticityPlaneStrain::cauchyStrain :
        NULL;
    kernels[1] = ProjectKernels("cauchy_strain", strainKernel);

    assert(integrator);
    integrator->setKernelsDerivedField(kernels);

    PYLITH_METHOD_END;
} // _setKernelsDerivedField


// End of file
