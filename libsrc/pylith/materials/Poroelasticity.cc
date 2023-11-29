// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/materials/Poroelasticity.hh" // implementation of object methods

#include "pylith/materials/RheologyPoroelasticity.hh" // HASA RheologyPoroelasticity
#include "pylith/materials/AuxiliaryFactoryPoroelastic.hh" // USES AuxiliaryFactory
#include "pylith/materials/DerivedFactoryPoroelasticity.hh" // USES DerivedFactoryPoroelasticity
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/fekernels/Poroelasticity.hh" // USES Poroelasticity kernels
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/DispVel.hh" // USES DispVel kernels

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;
typedef pylith::feassemble::Integrator::EquationPart EquationPart;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::Poroelasticity::Poroelasticity(void) :
    _useBodyForce(false),
    _useReferenceState(false),
    _useSourceDensity(false),
    _useStateVars(false),
    _rheology(NULL),
    _derivedFactory(new pylith::materials::DerivedFactoryPoroelasticity) {
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
// Update Fields?
void
pylith::materials::Poroelasticity::useStateVars(const bool value) {
    PYLITH_COMPONENT_DEBUG("useStateVars(value=" << value << ")");
    _useStateVars = value;
} // useStateVars


// ---------------------------------------------------------------------------------------------------------------------
// Update fields?
bool
pylith::materials::Poroelasticity::useStateVars(void) const {
    return _useStateVars;
} // useStateVars


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

    // Verify solution contains required fields.
    std::string reason = "material 'Poroelasticity'.";
    size_t numRequired = 0;
    const size_t maxRequired = 6;
    pylith::string_vector requiredFields(maxRequired);
    requiredFields[numRequired++] = "displacement";
    requiredFields[numRequired++] = "pressure";

    switch (_formulation) {
    case QUASISTATIC:
        requiredFields[numRequired++] = "trace_strain";
        if (_useStateVars) {
            reason = "material 'Poroelasticity' with state variables.";
            requiredFields[numRequired++] = "velocity";
            requiredFields[numRequired++] = "pressure_t";
            requiredFields[numRequired++] = "trace_strain_t";
        } // if
        break;
    case DYNAMIC:
    case DYNAMIC_IMEX:
        reason = "material 'Poroelasticity' with inertia.";
        requiredFields[numRequired++] = "velocity";
    } // switch
    requiredFields.resize(numRequired);

    pylith::topology::FieldOps::checkSubfieldsExist(requiredFields, reason, solution);

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::materials::Poroelasticity::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    pylith::feassemble::IntegratorDomain* integrator = new pylith::feassemble::IntegratorDomain(this);assert(integrator);
    integrator->setLabelName(getLabelName());
    integrator->setLabelValue(getLabelValue());

    _setKernelsResidual(integrator, solution);
    _setKernelsJacobian(integrator, solution);
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
    auxiliaryFactory->initialize(auxiliaryField, *_normalizer, domainMesh.getDimension());

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
    _derivedFactory->initialize(derivedField, *_normalizer, domainMesh.getDimension());
    _derivedFactory->addSubfields();

    derivedField->subfieldsSetup();
    derivedField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *derivedField);
    derivedField->allocate();
    derivedField->createOutputVector();

    PYLITH_METHOD_RETURN(derivedField);
} // createDerivedField


// ------------------------------------------------------------------------------------------------
// Get default PETSc solver options appropriate for material.
pylith::utils::PetscOptions*
pylith::materials::Poroelasticity::getSolverDefaults(const bool isParallel,
                                                     const bool hasFault) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getSolverDefaults(isParallel="<<isParallel<<", hasFault="<<hasFault<<")");

    pylith::utils::PetscOptions* options = new pylith::utils::PetscOptions();assert(options);

    switch (_formulation) {
    case pylith::problems::Physics::QUASISTATIC:
        options->add("-ts_type", "beuler");

        if (!hasFault) {
            if (!isParallel) {
                options->add("-pc_type", "lu");
            } else {
#if 1
                options->add("-pc_type", "ml");
#else
                options->add("-pc_type", "gamg");
                options->add("-mg_levels_pc_type", "sor");
                options->add("-mg_levels_ksp_type", "richardson");
#endif
            } // if/else
        } // if
        break;
    case pylith::problems::Physics::DYNAMIC:
        break;
    case pylith::problems::Physics::DYNAMIC_IMEX:
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation '" << _formulation << "'.");
    } // switch

    PYLITH_METHOD_RETURN(options);
} // getSolverDefaults


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
// Set kernels for residual.
void
pylith::materials::Poroelasticity::_setKernelsResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                       const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsResidual(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();

    const int bitBodyForce = _useBodyForce ? 0x1 : 0x0;
    const int bitGravity = _gravityField ? 0x2 : 0x0;
    const int bitSourceDensity = _useSourceDensity ? 0x4 : 0x0;
    const int bitUse = bitBodyForce | bitGravity | bitSourceDensity;

    PetscPointFunc r0 = NULL;
    switch (bitUse) {
    case 0x0:
        break;
    case 0x1:
        r0 = pylith::fekernels::Poroelasticity::g0v_bodyforce;
        break;
    case 0x2:
        r0 = pylith::fekernels::Poroelasticity::g0v_grav;
        break;
    case 0x4:
        break;
    case 0x3:
        r0 = pylith::fekernels::Poroelasticity::g0v_grav_bodyforce;
        break;
    case 0x5:
        r0 = pylith::fekernels::Poroelasticity::g0v_bodyforce;
        break;
    case 0x6:
        r0 = pylith::fekernels::Poroelasticity::g0v_grav;
        break;
    case 0x7:
        r0 = pylith::fekernels::Poroelasticity::g0v_grav_bodyforce;
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown case (bitUse=" << bitUse << ") for residual kernels.");
    } // switch

    std::vector<ResidualKernels> kernels;
    switch (_formulation) {
    case QUASISTATIC: {
        if (!_useStateVars) {
            // Displacement
            const PetscPointFunc f0u = r0;
            const PetscPointFunc f1u = _rheology->getKernelf1u_implicit(coordsys);

            // Pressure
            PetscPointFunc f0p = _rheology->getKernelf0p_implicit(coordsys, _useBodyForce, _gravityField, _useSourceDensity);
            PetscPointFunc f1p = _rheology->getKernelf1p_implicit(coordsys, _useBodyForce, _gravityField); // darcy
                                                                                                           // velocity

            // Volumetric Strain
            const PetscPointFunc f0e = pylith::fekernels::Poroelasticity::f0e;
            const PetscPointFunc f1e = NULL;

            kernels.resize(3);
            kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::LHS, f0u, f1u);
            kernels[1] = ResidualKernels("pressure", pylith::feassemble::Integrator::LHS, f0p, f1p);
            kernels[2] = ResidualKernels("trace_strain", pylith::feassemble::Integrator::LHS, f0e, f1e);
        } else {
            // Displacement
            PetscPointFunc f0u = r0;
            const PetscPointFunc f1u = _rheology->getKernelf1u_implicit(coordsys);

            // Pressure
            PetscPointFunc f0p = _rheology->getKernelf0p_implicit(coordsys, _useBodyForce, _gravityField, _useSourceDensity);
            PetscPointFunc f1p = _rheology->getKernelf1p_implicit(coordsys, _useBodyForce, _gravityField);

            // Volumetric Strain
            const PetscPointFunc f0e = pylith::fekernels::Poroelasticity::f0e;
            const PetscPointFunc f1e = NULL;

            // Velocity
            const PetscPointFunc f0v = pylith::fekernels::Poroelasticity::f0v_implicit;
            const PetscPointFunc f1v = NULL;

            // Time derivative of pressure
            const PetscPointFunc f0pdot = pylith::fekernels::Poroelasticity::f0pdot;
            const PetscPointFunc f1pdot = NULL;

            // Time derivative of Volumetric Strain
            const PetscPointFunc f0edot = pylith::fekernels::Poroelasticity::f0edot;
            const PetscPointFunc f1edot = NULL;

            kernels.resize(6);
            kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::LHS, f0u, f1u);
            kernels[1] = ResidualKernels("pressure", pylith::feassemble::Integrator::LHS, f0p, f1p);
            kernels[2] = ResidualKernels("trace_strain", pylith::feassemble::Integrator::LHS, f0e, f1e);
            kernels[3] = ResidualKernels("velocity", pylith::feassemble::Integrator::LHS, f0v, f1v);
            kernels[4] = ResidualKernels("pressure_t", pylith::feassemble::Integrator::LHS, f0pdot, f1pdot);
            kernels[5] = ResidualKernels("trace_strain_t", pylith::feassemble::Integrator::LHS, f0edot, f1edot);
        }
        break;
    } // QUASISTATIC
    case DYNAMIC_IMEX:
    case DYNAMIC: {
        // Displacement
        const PetscPointFunc f0u = pylith::fekernels::DispVel::f0u;
        const PetscPointFunc f1u = NULL;
        const PetscPointFunc g0u = pylith::fekernels::Poroelasticity::g0u;
        const PetscPointFunc g1u = NULL;

        // Pressure
        const PetscPointFunc f0p = _rheology->getKernelf0p_explicit(coordsys);
        const PetscPointFunc f1p = NULL;
        const PetscPointFunc g0p = _rheology->getKernelg0p(coordsys, _useBodyForce, _gravityField, _useSourceDensity);
        const PetscPointFunc g1p = _rheology->getKernelg1p_explicit(coordsys, _gravityField); // Darcy velocity

        // Velocity
        const PetscPointFunc f0v = pylith::fekernels::Poroelasticity::f0v_explicit;
        const PetscPointFunc f1v = NULL;
        const PetscPointFunc g0v = r0;
        const PetscPointFunc g1v = _rheology->getKernelg1v_explicit(coordsys);

        kernels.resize(6);
        kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::LHS, f0u, f1u);
        kernels[1] = ResidualKernels("displacement", pylith::feassemble::Integrator::RHS, g0u, g1u);
        kernels[2] = ResidualKernels("pressure", pylith::feassemble::Integrator::LHS, f0p, f1p);
        kernels[3] = ResidualKernels("pressure", pylith::feassemble::Integrator::RHS, g0p, g1p);
        kernels[4] = ResidualKernels("velocity", pylith::feassemble::Integrator::LHS, f0v, f1v);
        kernels[5] = ResidualKernels("velocity", pylith::feassemble::Integrator::RHS, g0v, g1v);
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    // Add any MMS body force kernels.
    kernels.insert(kernels.end(), _mmsBodyForceKernels.begin(), _mmsBodyForceKernels.end());

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsResidual


// ----------------------------------------------------------------------
// Set kernels for Jacobian.
void
pylith::materials::Poroelasticity::_setKernelsJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                       const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsJacobian(integrator="<<integrator<<",solution="<<solution.getLabel()<<")");
    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    std::vector<JacobianKernels> kernels(9);

    switch (_formulation) {
    case QUASISTATIC: {
        const EquationPart equationPart = pylith::feassemble::Integrator::LHS;
        if (!_useStateVars) {
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

            kernels[0] = JacobianKernels("displacement", "displacement", equationPart, Jf0uu, Jf1uu, Jf2uu, Jf3uu);
            kernels[1] = JacobianKernels("displacement", "pressure",     equationPart, Jf0up, Jf1up, Jf2up, Jf3up);
            kernels[2] = JacobianKernels("displacement", "trace_strain", equationPart, Jf0ue, Jf1ue, Jf2ue, Jf3ue);
            kernels[3] = JacobianKernels("pressure",     "displacement", equationPart, Jf0pu, Jf1pu, Jf2pu, Jf3pu);
            kernels[4] = JacobianKernels("pressure",     "pressure",     equationPart, Jf0pp, Jf1pp, Jf2pp, Jf3pp);
            kernels[5] = JacobianKernels("pressure",     "trace_strain", equationPart, Jf0pe, Jf1pe, Jf2pe, Jf3pe);
            kernels[6] = JacobianKernels("trace_strain", "displacement", equationPart, Jf0eu, Jf1eu, Jf2eu, Jf3eu);
            kernels[7] = JacobianKernels("trace_strain", "pressure",     equationPart, Jf0ep, Jf1ep, Jf2ep, Jf3ep);
            kernels[8] = JacobianKernels("trace_strain", "trace_strain", equationPart, Jf0ee, Jf1ee, Jf2ee, Jf3ee);
        } else {
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

            const PetscPointJac Jf0uv = NULL;
            const PetscPointJac Jf1uv = NULL;
            const PetscPointJac Jf2uv = NULL;
            const PetscPointJac Jf3uv = NULL;

            const PetscPointJac Jf0updot = NULL;
            const PetscPointJac Jf1updot = NULL;
            const PetscPointJac Jf2updot = NULL;
            const PetscPointJac Jf3updot = NULL;

            const PetscPointJac Jf0uedot = NULL;
            const PetscPointJac Jf1uedot = NULL;
            const PetscPointJac Jf2uedot = NULL;
            const PetscPointJac Jf3uedot = NULL;

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

            const PetscPointJac Jf0pv = NULL;
            const PetscPointJac Jf1pv = NULL;
            const PetscPointJac Jf2pv = NULL;
            const PetscPointJac Jf3pv = NULL;

            const PetscPointJac Jf0ppdot = _rheology->getKernelJf0ppdot(coordsys);
            const PetscPointJac Jf1ppdot = NULL;
            const PetscPointJac Jf2ppdot = NULL;
            const PetscPointJac Jf3ppdot = NULL;

            const PetscPointJac Jf0pedot = _rheology->getKernelJf0pedot(coordsys);
            const PetscPointJac Jf1pedot = NULL;
            const PetscPointJac Jf2pedot = NULL;
            const PetscPointJac Jf3pedot = NULL;

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

            const PetscPointJac Jf0ev = NULL;
            const PetscPointJac Jf1ev = NULL;
            const PetscPointJac Jf2ev = NULL;
            const PetscPointJac Jf3ev = NULL;

            const PetscPointJac Jf0epdot = NULL;
            const PetscPointJac Jf1epdot = NULL;
            const PetscPointJac Jf2epdot = NULL;
            const PetscPointJac Jf3epdot = NULL;

            const PetscPointJac Jf0eedot = NULL;
            const PetscPointJac Jf1eedot = NULL;
            const PetscPointJac Jf2eedot = NULL;
            const PetscPointJac Jf3eedot = NULL;

            const PetscPointJac Jf0vu = pylith::fekernels::Poroelasticity::Jf0vu_implicit;
            const PetscPointJac Jf1vu = NULL;
            const PetscPointJac Jf2vu = NULL;
            const PetscPointJac Jf3vu = NULL;

            const PetscPointJac Jf0vp = NULL;
            const PetscPointJac Jf1vp = NULL;
            const PetscPointJac Jf2vp = NULL;
            const PetscPointJac Jf3vp = NULL;

            const PetscPointJac Jf0ve = NULL;
            const PetscPointJac Jf1ve = NULL;
            const PetscPointJac Jf2ve = NULL;
            const PetscPointJac Jf3ve = NULL;

            const PetscPointJac Jf0vv = pylith::fekernels::Poroelasticity::Jf0vv_implicit;
            const PetscPointJac Jf1vv = NULL;
            const PetscPointJac Jf2vv = NULL;
            const PetscPointJac Jf3vv = NULL;

            const PetscPointJac Jf0vpdot = NULL;
            const PetscPointJac Jf1vpdot = NULL;
            const PetscPointJac Jf2vpdot = NULL;
            const PetscPointJac Jf3vpdot = NULL;

            const PetscPointJac Jf0vedot = NULL;
            const PetscPointJac Jf1vedot = NULL;
            const PetscPointJac Jf2vedot = NULL;
            const PetscPointJac Jf3vedot = NULL;

            const PetscPointJac Jf0pdotu = NULL;
            const PetscPointJac Jf1pdotu = NULL;
            const PetscPointJac Jf2pdotu = NULL;
            const PetscPointJac Jf3pdotu = NULL;

            const PetscPointJac Jf0pdotp = pylith::fekernels::Poroelasticity::Jf0pdotp;
            const PetscPointJac Jf1pdotp = NULL;
            const PetscPointJac Jf2pdotp = NULL;
            const PetscPointJac Jf3pdotp = NULL;

            const PetscPointJac Jf0pdote = NULL;
            const PetscPointJac Jf1pdote = NULL;
            const PetscPointJac Jf2pdote = NULL;
            const PetscPointJac Jf3pdote = NULL;

            const PetscPointJac Jf0pdotv = NULL;
            const PetscPointJac Jf1pdotv = NULL;
            const PetscPointJac Jf2pdotv = NULL;
            const PetscPointJac Jf3pdotv = NULL;

            const PetscPointJac Jf0pdotpdot = pylith::fekernels::Poroelasticity::Jf0pdotpdot;
            const PetscPointJac Jf1pdotpdot = NULL;
            const PetscPointJac Jf2pdotpdot = NULL;
            const PetscPointJac Jf3pdotpdot = NULL;

            const PetscPointJac Jf0pdotedot = NULL;
            const PetscPointJac Jf1pdotedot = NULL;
            const PetscPointJac Jf2pdotedot = NULL;
            const PetscPointJac Jf3pdotedot = NULL;

            const PetscPointJac Jf0edotu = NULL;
            const PetscPointJac Jf1edotu = NULL;
            const PetscPointJac Jf2edotu = NULL;
            const PetscPointJac Jf3edotu = NULL;

            const PetscPointJac Jf0edotp = NULL;
            const PetscPointJac Jf1edotp = NULL;
            const PetscPointJac Jf2edotp = NULL;
            const PetscPointJac Jf3edotp = NULL;

            const PetscPointJac Jf0edote = pylith::fekernels::Poroelasticity::Jf0edote;
            const PetscPointJac Jf1edote = NULL;
            const PetscPointJac Jf2edote = NULL;
            const PetscPointJac Jf3edote = NULL;

            const PetscPointJac Jf0edotv = NULL;
            const PetscPointJac Jf1edotv = NULL;
            const PetscPointJac Jf2edotv = NULL;
            const PetscPointJac Jf3edotv = NULL;

            const PetscPointJac Jf0edotpdot = NULL;
            const PetscPointJac Jf1edotpdot = NULL;
            const PetscPointJac Jf2edotpdot = NULL;
            const PetscPointJac Jf3edotpdot = NULL;

            const PetscPointJac Jf0edotedot = pylith::fekernels::Poroelasticity::Jf0edotedot;
            const PetscPointJac Jf1edotedot = NULL;
            const PetscPointJac Jf2edotedot = NULL;
            const PetscPointJac Jf3edotedot = NULL;

            kernels.resize(36);
            kernels[0] = JacobianKernels("displacement", "displacement", equationPart, Jf0uu, Jf1uu, Jf2uu, Jf3uu);
            kernels[1] = JacobianKernels("displacement", "pressure", equationPart, Jf0up, Jf1up, Jf2up, Jf3up);
            kernels[2] = JacobianKernels("displacement", "trace_strain", equationPart, Jf0ue, Jf1ue, Jf2ue, Jf3ue);
            kernels[3] = JacobianKernels("displacement", "velocity", equationPart, Jf0uv, Jf1uv, Jf2uv, Jf3uv);
            kernels[4] = JacobianKernels("displacement", "pressure_t", equationPart, Jf0updot, Jf1updot, Jf2updot, Jf3updot);
            kernels[5] = JacobianKernels("displacement", "trace_strain_t", equationPart, Jf0uedot, Jf1uedot, Jf2uedot, Jf3uedot);
            kernels[6] = JacobianKernels("pressure", "displacement", equationPart, Jf0pu, Jf1pu, Jf2pu, Jf3pu);
            kernels[7] = JacobianKernels("pressure", "pressure", equationPart, Jf0pp, Jf1pp, Jf2pp, Jf3pp);
            kernels[8] = JacobianKernels("pressure", "trace_strain", equationPart, Jf0pe, Jf1pe, Jf2pe, Jf3pe);
            kernels[9] = JacobianKernels("pressure", "velocity", equationPart, Jf0pv, Jf1pv, Jf2pv, Jf3pv);
            kernels[10] = JacobianKernels("pressure", "pressure_t", equationPart, Jf0ppdot, Jf1ppdot, Jf2ppdot, Jf3ppdot);
            kernels[11] = JacobianKernels("pressure", "trace_strain_t", equationPart, Jf0pedot, Jf1pedot, Jf2pedot, Jf3pedot);
            kernels[12] = JacobianKernels("trace_strain", "displacement", equationPart, Jf0eu, Jf1eu, Jf2eu, Jf3eu);
            kernels[13] = JacobianKernels("trace_strain", "pressure", equationPart, Jf0ep, Jf1ep, Jf2ep, Jf3ep);
            kernels[14] = JacobianKernels("trace_strain", "trace_strain", equationPart, Jf0ee, Jf1ee, Jf2ee, Jf3ee);
            kernels[15] = JacobianKernels("trace_strain", "velocity", equationPart, Jf0ev, Jf1ev, Jf2ev, Jf3ev);
            kernels[16] = JacobianKernels("trace_strain", "pressure_t", equationPart, Jf0epdot, Jf1epdot, Jf2epdot, Jf3epdot);
            kernels[17] = JacobianKernels("trace_strain", "trace_strain_t", equationPart, Jf0eedot, Jf1eedot, Jf2eedot, Jf3eedot);
            kernels[18] = JacobianKernels("velocity", "displacement", equationPart, Jf0vu, Jf1vu, Jf2vu, Jf3vu);
            kernels[19] = JacobianKernels("velocity", "pressure", equationPart, Jf0vp, Jf1vp, Jf2vp, Jf3vp);
            kernels[20] = JacobianKernels("velocity", "trace_strain", equationPart, Jf0ve, Jf1ve, Jf2ve, Jf3ve);
            kernels[21] = JacobianKernels("velocity", "velocity", equationPart, Jf0vv, Jf1vv, Jf2vv, Jf3vv);
            kernels[22] = JacobianKernels("velocity", "pressure_t", equationPart, Jf0vpdot, Jf1vpdot, Jf2vpdot, Jf3vpdot);
            kernels[23] = JacobianKernels("velocity", "trace_strain_t", equationPart, Jf0vedot, Jf1vedot, Jf2vedot, Jf3vedot);
            kernels[24] = JacobianKernels("pressure_t", "displacement", equationPart, Jf0pdotu, Jf1pdotu, Jf2pdotu, Jf3pdotu);
            kernels[25] = JacobianKernels("pressure_t", "pressure", equationPart, Jf0pdotp, Jf1pdotp, Jf2pdotp, Jf3pdotp);
            kernels[26] = JacobianKernels("pressure_t", "trace_strain", equationPart, Jf0pdote, Jf1pdote, Jf2pdote, Jf3pdote);
            kernels[27] = JacobianKernels("pressure_t", "velocity", equationPart, Jf0pdotv, Jf1pdotv, Jf2pdotv, Jf3pdotv);
            kernels[28] = JacobianKernels("pressure_t", "pressure_t", equationPart, Jf0pdotpdot, Jf1pdotpdot, Jf2pdotpdot, Jf3pdotpdot);
            kernels[29] = JacobianKernels("pressure_t", "trace_strain_t", equationPart, Jf0pdotedot, Jf1pdotedot, Jf2pdotedot, Jf3pdotedot);
            kernels[30] = JacobianKernels("trace_strain_t", "displacement", equationPart, Jf0edotu, Jf1edotu, Jf2edotu, Jf3edotu);
            kernels[31] = JacobianKernels("trace_strain_t", "pressure", equationPart, Jf0edotp, Jf1edotp, Jf2edotp, Jf3edotp);
            kernels[32] = JacobianKernels("trace_strain_t", "trace_strain", equationPart, Jf0edote, Jf1edote, Jf2edote, Jf3edote);
            kernels[33] = JacobianKernels("trace_strain_t", "velocity", equationPart, Jf0edotv, Jf1edotv, Jf2edotv, Jf3edotv);
            kernels[34] = JacobianKernels("trace_strain_t", "pressure_t", equationPart, Jf0edotpdot, Jf1edotpdot, Jf2edotpdot, Jf3edotpdot);
            kernels[35] = JacobianKernels("trace_strain_t", "trace_strain_t", equationPart, Jf0edotedot, Jf1edotedot, Jf2edotedot, Jf3edotedot);
        }
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

        const PetscPointJac Jf0vv = pylith::fekernels::Poroelasticity::Jf0vv_explicit;
        const PetscPointJac Jf1vv = NULL;
        const PetscPointJac Jf2vv = NULL;
        const PetscPointJac Jf3vv = NULL;

        integrator->setLHSJacobianTriggers(pylith::feassemble::Integrator::NEW_JACOBIAN_TIME_STEP_CHANGE);

        const EquationPart equationPart = pylith::feassemble::Integrator::LHS_LUMPED_INV;
        kernels[0] = JacobianKernels("displacement",  "displacement", equationPart, Jf0uu, Jf1uu, Jf2uu, Jf3uu);
        kernels[1] = JacobianKernels("displacement",  "pressure",     equationPart, Jf0up, Jf1up, Jf2up, Jf3up);
        kernels[2] = JacobianKernels("displacement",  "velocity",     equationPart, Jf0uv, Jf1uv, Jf2uv, Jf3uv);
        kernels[3] = JacobianKernels("pressure",      "displacement", equationPart, Jf0pu, Jf1pu, Jf2pu, Jf3pu);
        kernels[4] = JacobianKernels("pressure",      "pressure",     equationPart, Jf0pp, Jf1pp, Jf2pp, Jf3pp);
        kernels[5] = JacobianKernels("pressure",      "velocity",     equationPart, Jf0pv, Jf1pv, Jf2pv, Jf3pv);
        kernels[6] = JacobianKernels("velocity",      "displacement", equationPart, Jf0vu, Jf1vu, Jf2vu, Jf3vu);
        kernels[7] = JacobianKernels("velocity",      "pressure",     equationPart, Jf0vp, Jf1vp, Jf2vp, Jf3vp);
        kernels[8] = JacobianKernels("velocity",      "velocity",     equationPart, Jf0vv, Jf1vv, Jf2vv, Jf3vv);
        break;
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    assert(integrator);
    integrator->setKernelsJacobian(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for computing updated state variables in auxiliary field.
void
pylith::materials::Poroelasticity::_setKernelsUpdateStateVars(pylith::feassemble::IntegratorDomain* integrator,
                                                              const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsUpdateStateVars(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    std::vector<ProjectKernels> kernels;

    switch (_formulation) {
    case QUASISTATIC: {
        _rheology->addKernelsUpdateStateVarsImplicit(&kernels, coordsys, _useStateVars);

        break;
    } // QUASISTATIC
    case DYNAMIC_IMEX:
    case DYNAMIC: {
        _rheology->addKernelsUpdateStateVarsExplicit(&kernels, coordsys, _useStateVars);

        break;
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

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

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    std::vector<ProjectKernels> kernels(4);
    kernels[0] = ProjectKernels("cauchy_stress", _rheology->getKernelCauchyStressVector(coordsys));

    const int spaceDim = coordsys->getSpaceDim();
    const PetscPointFunc strainKernel =
        (3 == spaceDim) ? pylith::fekernels::Elasticity3D::infinitesimalStrain_asVector :
        (2 == spaceDim) ? pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain_asVector :
        NULL;
    kernels[1] = ProjectKernels("cauchy_strain", strainKernel);

    const PetscPointFunc bulkDensity = pylith::fekernels::Poroelasticity::bulkDensity_asScalar;

    kernels[2] = ProjectKernels("bulk_density", bulkDensity);

    kernels[3] = ProjectKernels("water_content", _rheology->getKernelWaterContent(coordsys));

    assert(integrator);
    integrator->setKernelsDerivedField(kernels);

    PYLITH_METHOD_END;
} // _setKernelsDerivedField


// End of file
