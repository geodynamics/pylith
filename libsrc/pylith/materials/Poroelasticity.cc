// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

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
#include "pylith/scales/Scales.hh" // USES Scales

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
    integrator->createLabelDS(solution, solution.getMesh().getDimension());

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
    auxiliaryField->setLabel("auxiliary field");

    assert(_rheology);
    pylith::materials::AuxiliaryFactoryPoroelasticity* auxiliaryFactory = _rheology->getAuxiliaryFactory();assert(auxiliaryFactory);

    assert(_scales);
    auxiliaryFactory->initialize(auxiliaryField, *_scales, domainMesh.getDimension());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    // ---------------------------------
    // Required Auxiliary
    auxiliaryFactory->addSolidDensity(); // 0 Rock Density
    auxiliaryFactory->addFluidDensity(); // 1 Fluid Density
    auxiliaryFactory->addFluidViscosity(); // 2 Fluid Viscosity
    auxiliaryFactory->addPorosity(); // 3 Porosity

    // ---------------------------------
    // Optional Auxiliary
    if (_useBodyForce) {
        auxiliaryFactory->addBodyForce();
    } // if
    if (_gravityField) {
        auxiliaryFactory->addGravityField(_gravityField);
    } // if
    if (_useSourceDensity) {
        auxiliaryFactory->addSourceDensity();
    } // if
    _rheology->addAuxiliarySubfields();

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    assert(auxiliaryFactory);
    auxiliaryFactory->setValuesFromDB();

    pythia::journal::debug_t debug("poroelasticity.view_auxiliary_field");
    if (debug.state()) {
        auxiliaryField->view("Poroelasticity auxiliary field");
    } // if

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
    derivedField->setLabel("derived field");

    assert(_scales);
    _derivedFactory->initialize(derivedField, *_scales, domainMesh.getDimension());
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
        options->add("-ts_exact_final_time", "matchstep");

        options->add("-pc_type", "fieldsplit");
        options->add("-pc_fieldsplit_type", "multiplicative");
        options->add("-pc_fieldsplit_0_fields", "2");
        options->add("-pc_fieldsplit_1_fields", "1");
        options->add("-pc_fieldsplit_2_fields", "0");
        options->add("-fieldsplit_trace_strain_pc_type", "bjacobi");
        options->add("-fieldsplit_pressure_pc_type", "bjacobi");

        if (!isParallel) {
            options->add("-fieldsplit_displacement_pc_type", "lu");
        } else {
            options->add("-fieldsplit_displacement_pc_type", "ml");
            options->add("-fieldsplit_displacement_ksp_type", "gmres");
        } // if/else

        if (_useStateVars) {
            options->add("-fieldsplit_velocity_pc_type", "bjacobi");
            options->add("-fieldsplit_pressure_t_pc_type", "bjacobi");
            options->add("-fieldsplit_trace_strain_t_pc_type", "bjacobi");
        } // if

        if (!hasFault) {
            if (_useStateVars) {
                options->add("-pc_fieldsplit_3_fields", "3");
                options->add("-pc_fieldsplit_4_fields", "4");
                options->add("-pc_fieldsplit_5_fields", "5");
            } // if

        } else {
            options->add("-pc_fieldsplit_3_fields", "3");
            options->add("-fieldsplit_lagrange_multiplier_fault_pc_type", "jacobi");

            options->add("-dm_reorder_section");
            options->add("-dm_reorder_section_type", "cohesive");
            options->add("-fieldsplit_displacement_pc_type", "ml");
            options->add("-fieldsplit_displacement_ksp_type", "gmres");
            options->add("-fieldsplit_displacement_mg_fine_pc_type", "vpbjacobi");

            if (_useStateVars) {
                options->add("-pc_fieldsplit_4_fields", "4");
                options->add("-pc_fieldsplit_5_fields", "5");
                options->add("-pc_fieldsplit_6_fields", "6");
            } // if
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
    const int bitUse = bitBodyForce | bitGravity;

    PetscPointFn* r0 = NULL;
    switch (bitUse) {
    case 0x0:
        break;
    case 0x1:
        r0 = pylith::fekernels::Poroelasticity::g0v_bodyforce;
        break;
    case 0x2:
        r0 = pylith::fekernels::Poroelasticity::g0v_grav;
        break;
    case 0x3:
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
            const PetscPointFn* f0u = r0;
            const PetscPointFn* f1u = _rheology->getKernelf1u_implicit(coordsys);

            // Pressure
            PetscPointFn* f0p = _rheology->getKernelf0p_implicit(coordsys, _useBodyForce, _gravityField, _useSourceDensity);
            PetscPointFn* f1p = _rheology->getKernelf1p_implicit(coordsys, _useBodyForce, _gravityField); // darcy
                                                                                                          // velocity

            // Volumetric Strain
            const PetscPointFn* f0e = pylith::fekernels::Poroelasticity::f0e;
            const PetscPointFn* f1e = NULL;

            kernels.resize(3);
            kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::LHS, f0u, f1u);
            kernels[1] = ResidualKernels("pressure", pylith::feassemble::Integrator::LHS, f0p, f1p);
            kernels[2] = ResidualKernels("trace_strain", pylith::feassemble::Integrator::LHS, f0e, f1e);
        } else {
            // Displacement
            PetscPointFn* f0u = r0;
            const PetscPointFn* f1u = _rheology->getKernelf1u_implicit(coordsys);

            // Pressure
            PetscPointFn* f0p = _rheology->getKernelf0p_implicit(coordsys, _useBodyForce, _gravityField, _useSourceDensity);
            PetscPointFn* f1p = _rheology->getKernelf1p_implicit(coordsys, _useBodyForce, _gravityField);

            // Volumetric Strain
            const PetscPointFn* f0e = pylith::fekernels::Poroelasticity::f0e;
            const PetscPointFn* f1e = NULL;

            // Velocity
            const PetscPointFn* f0v = pylith::fekernels::Poroelasticity::f0v_implicit;
            const PetscPointFn* f1v = NULL;

            // Time derivative of pressure
            const PetscPointFn* f0pdot = pylith::fekernels::Poroelasticity::f0pdot;
            const PetscPointFn* f1pdot = NULL;

            // Time derivative of Volumetric Strain
            const PetscPointFn* f0edot = pylith::fekernels::Poroelasticity::f0edot;
            const PetscPointFn* f1edot = NULL;

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
        const PetscPointFn* f0u = pylith::fekernels::DispVel::f0u;
        const PetscPointFn* f1u = NULL;
        const PetscPointFn* g0u = pylith::fekernels::Poroelasticity::g0u;
        const PetscPointFn* g1u = NULL;

        // Pressure
        const PetscPointFn* f0p = _rheology->getKernelf0p_explicit(coordsys);
        const PetscPointFn* f1p = NULL;
        const PetscPointFn* g0p = _rheology->getKernelg0p(coordsys, _useBodyForce, _gravityField, _useSourceDensity);
        const PetscPointFn* g1p = _rheology->getKernelg1p_explicit(coordsys, _gravityField); // Darcy velocity

        // Velocity
        const PetscPointFn* f0v = pylith::fekernels::Poroelasticity::f0v_explicit;
        const PetscPointFn* f1v = NULL;
        const PetscPointFn* g0v = r0;
        const PetscPointFn* g1v = _rheology->getKernelg1v_explicit(coordsys);

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
    std::vector<JacobianKernels> kernels(7);

    integrator->setLHSJacobianTriggers(pylith::feassemble::Integrator::NEW_JACOBIAN_TIME_STEP_CHANGE);

    switch (_formulation) {
    case QUASISTATIC: {
        const EquationPart equationPart = pylith::feassemble::Integrator::LHS;
        if (!_useStateVars) {
            PetscPointJacFn* Jf0uu = NULL;
            PetscPointJacFn* Jf1uu = NULL;
            PetscPointJacFn* Jf2uu = NULL;
            PetscPointJacFn* Jf3uu = _rheology->getKernelJf3uu(coordsys);

            PetscPointJacFn* Jf0up = NULL;
            PetscPointJacFn* Jf1up = NULL;
            PetscPointJacFn* Jf2up = _rheology->getKernelJf2up(coordsys);
            PetscPointJacFn* Jf3up = NULL;

            PetscPointJacFn* Jf0ue = NULL;
            PetscPointJacFn* Jf1ue = NULL;
            PetscPointJacFn* Jf2ue = _rheology->getKernelJf2ue(coordsys);
            PetscPointJacFn* Jf3ue = NULL;

            PetscPointJacFn* Jf0pp = _rheology->getKernelJf0pp(coordsys);
            PetscPointJacFn* Jf1pp = NULL;
            PetscPointJacFn* Jf2pp = NULL;
            PetscPointJacFn* Jf3pp = _rheology->getKernelJf3pp(coordsys);

            PetscPointJacFn* Jf0pe = _rheology->getKernelJf0pe(coordsys);
            PetscPointJacFn* Jf1pe = NULL;
            PetscPointJacFn* Jf2pe = NULL;
            PetscPointJacFn* Jf3pe = NULL;

            PetscPointJacFn* Jf0eu = NULL;
            PetscPointJacFn* Jf1eu = pylith::fekernels::Poroelasticity::Jf1eu;
            PetscPointJacFn* Jf2eu = NULL;
            PetscPointJacFn* Jf3eu = NULL;

            PetscPointJacFn* Jf0ee = pylith::fekernels::Poroelasticity::Jf0ee;
            PetscPointJacFn* Jf1ee = NULL;
            PetscPointJacFn* Jf2ee = NULL;
            PetscPointJacFn* Jf3ee = NULL;

            kernels[0] = JacobianKernels("displacement", "displacement", equationPart, Jf0uu, Jf1uu, Jf2uu, Jf3uu);
            kernels[1] = JacobianKernels("displacement", "pressure",     equationPart, Jf0up, Jf1up, Jf2up, Jf3up);
            kernels[2] = JacobianKernels("displacement", "trace_strain", equationPart, Jf0ue, Jf1ue, Jf2ue, Jf3ue);
            kernels[3] = JacobianKernels("pressure",     "pressure",     equationPart, Jf0pp, Jf1pp, Jf2pp, Jf3pp);
            kernels[4] = JacobianKernels("pressure",     "trace_strain", equationPart, Jf0pe, Jf1pe, Jf2pe, Jf3pe);
            kernels[5] = JacobianKernels("trace_strain", "displacement", equationPart, Jf0eu, Jf1eu, Jf2eu, Jf3eu);
            kernels[6] = JacobianKernels("trace_strain", "trace_strain", equationPart, Jf0ee, Jf1ee, Jf2ee, Jf3ee);
        } else {
            PetscPointJacFn* Jf0uu = NULL;
            PetscPointJacFn* Jf1uu = NULL;
            PetscPointJacFn* Jf2uu = NULL;
            PetscPointJacFn* Jf3uu = _rheology->getKernelJf3uu(coordsys);

            PetscPointJacFn* Jf0up = NULL;
            PetscPointJacFn* Jf1up = NULL;
            PetscPointJacFn* Jf2up = _rheology->getKernelJf2up(coordsys);
            PetscPointJacFn* Jf3up = NULL;

            PetscPointJacFn* Jf0ue = NULL;
            PetscPointJacFn* Jf1ue = NULL;
            PetscPointJacFn* Jf2ue = _rheology->getKernelJf2ue(coordsys);
            PetscPointJacFn* Jf3ue = NULL;

            PetscPointJacFn* Jf0pp = _rheology->getKernelJf0pp(coordsys);
            PetscPointJacFn* Jf1pp = NULL;
            PetscPointJacFn* Jf2pp = NULL;
            PetscPointJacFn* Jf3pp = _rheology->getKernelJf3pp(coordsys);

            PetscPointJacFn* Jf0pe = _rheology->getKernelJf0pe(coordsys);
            PetscPointJacFn* Jf1pe = NULL;
            PetscPointJacFn* Jf2pe = NULL;
            PetscPointJacFn* Jf3pe = NULL;

            PetscPointJacFn* Jf0eu = NULL;
            PetscPointJacFn* Jf1eu = pylith::fekernels::Poroelasticity::Jf1eu;
            PetscPointJacFn* Jf2eu = NULL;
            PetscPointJacFn* Jf3eu = NULL;

            PetscPointJacFn* Jf0ee = pylith::fekernels::Poroelasticity::Jf0ee;
            PetscPointJacFn* Jf1ee = NULL;
            PetscPointJacFn* Jf2ee = NULL;
            PetscPointJacFn* Jf3ee = NULL;

            PetscPointJacFn* Jf0vu = pylith::fekernels::Poroelasticity::Jf0vu_implicit;
            PetscPointJacFn* Jf1vu = NULL;
            PetscPointJacFn* Jf2vu = NULL;
            PetscPointJacFn* Jf3vu = NULL;

            PetscPointJacFn* Jf0vv = pylith::fekernels::Poroelasticity::Jf0vv_implicit;
            PetscPointJacFn* Jf1vv = NULL;
            PetscPointJacFn* Jf2vv = NULL;
            PetscPointJacFn* Jf3vv = NULL;

            PetscPointJacFn* Jf0pdotp = pylith::fekernels::Poroelasticity::Jf0pdotp;
            PetscPointJacFn* Jf1pdotp = NULL;
            PetscPointJacFn* Jf2pdotp = NULL;
            PetscPointJacFn* Jf3pdotp = NULL;

            PetscPointJacFn* Jf0pdotpdot = pylith::fekernels::Poroelasticity::Jf0pdotpdot;
            PetscPointJacFn* Jf1pdotpdot = NULL;
            PetscPointJacFn* Jf2pdotpdot = NULL;
            PetscPointJacFn* Jf3pdotpdot = NULL;

            PetscPointJacFn* Jf0edote = pylith::fekernels::Poroelasticity::Jf0edote;
            PetscPointJacFn* Jf1edote = NULL;
            PetscPointJacFn* Jf2edote = NULL;
            PetscPointJacFn* Jf3edote = NULL;

            PetscPointJacFn* Jf0edotedot = pylith::fekernels::Poroelasticity::Jf0edotedot;
            PetscPointJacFn* Jf1edotedot = NULL;
            PetscPointJacFn* Jf2edotedot = NULL;
            PetscPointJacFn* Jf3edotedot = NULL;

            kernels.resize(13);
            kernels[0] = JacobianKernels("displacement", "displacement", equationPart, Jf0uu, Jf1uu, Jf2uu, Jf3uu);
            kernels[1] = JacobianKernels("displacement", "pressure", equationPart, Jf0up, Jf1up, Jf2up, Jf3up);
            kernels[2] = JacobianKernels("displacement", "trace_strain", equationPart, Jf0ue, Jf1ue, Jf2ue, Jf3ue);
            kernels[3] = JacobianKernels("pressure", "pressure", equationPart, Jf0pp, Jf1pp, Jf2pp, Jf3pp);
            kernels[4] = JacobianKernels("pressure", "trace_strain", equationPart, Jf0pe, Jf1pe, Jf2pe, Jf3pe);
            kernels[5] = JacobianKernels("trace_strain", "displacement", equationPart, Jf0eu, Jf1eu, Jf2eu, Jf3eu);
            kernels[6] = JacobianKernels("trace_strain", "trace_strain", equationPart, Jf0ee, Jf1ee, Jf2ee, Jf3ee);
            kernels[7] = JacobianKernels("velocity", "displacement", equationPart, Jf0vu, Jf1vu, Jf2vu, Jf3vu);
            kernels[8] = JacobianKernels("velocity", "velocity", equationPart, Jf0vv, Jf1vv, Jf2vv, Jf3vv);
            kernels[9] = JacobianKernels("pressure_t", "pressure", equationPart, Jf0pdotp, Jf1pdotp, Jf2pdotp, Jf3pdotp);
            kernels[10] = JacobianKernels("pressure_t", "pressure_t", equationPart, Jf0pdotpdot, Jf1pdotpdot, Jf2pdotpdot, Jf3pdotpdot);
            kernels[11] = JacobianKernels("trace_strain_t", "trace_strain", equationPart, Jf0edote, Jf1edote, Jf2edote, Jf3edote);
            kernels[12] = JacobianKernels("trace_strain_t", "trace_strain_t", equationPart, Jf0edotedot, Jf1edotedot, Jf2edotedot, Jf3edotedot);
        }
        break;
    } // QUASISTATIC
    case DYNAMIC_IMEX:
    case DYNAMIC: {
        PetscPointJacFn* Jf0uu = pylith::fekernels::DispVel::Jf0uu_stshift;
        PetscPointJacFn* Jf1uu = NULL;
        PetscPointJacFn* Jf2uu = NULL;
        PetscPointJacFn* Jf3uu = NULL;

        PetscPointJacFn* Jf0up = NULL;
        PetscPointJacFn* Jf1up = NULL;
        PetscPointJacFn* Jf2up = NULL;
        PetscPointJacFn* Jf3up = NULL;

        PetscPointJacFn* Jf0uv = NULL;
        PetscPointJacFn* Jf1uv = NULL;
        PetscPointJacFn* Jf2uv = NULL;
        PetscPointJacFn* Jf3uv = NULL;

        PetscPointJacFn* Jf0pu = NULL;
        PetscPointJacFn* Jf1pu = NULL;
        PetscPointJacFn* Jf2pu = NULL;
        PetscPointJacFn* Jf3pu = NULL;

        PetscPointJacFn* Jf0pp = _rheology->getKernelJf0pp(coordsys);
        PetscPointJacFn* Jf1pp = NULL;
        PetscPointJacFn* Jf2pp = NULL;
        PetscPointJacFn* Jf3pp = NULL;

        PetscPointJacFn* Jf0pv = NULL;
        PetscPointJacFn* Jf1pv = NULL;
        PetscPointJacFn* Jf2pv = NULL;
        PetscPointJacFn* Jf3pv = NULL;

        PetscPointJacFn* Jf0vu = NULL;
        PetscPointJacFn* Jf1vu = NULL;
        PetscPointJacFn* Jf2vu = NULL;
        PetscPointJacFn* Jf3vu = NULL;

        PetscPointJacFn* Jf0vp = NULL;
        PetscPointJacFn* Jf1vp = NULL;
        PetscPointJacFn* Jf2vp = NULL;
        PetscPointJacFn* Jf3vp = NULL;

        PetscPointJacFn* Jf0vv = pylith::fekernels::Poroelasticity::Jf0vv_explicit;
        PetscPointJacFn* Jf1vv = NULL;
        PetscPointJacFn* Jf2vv = NULL;
        PetscPointJacFn* Jf3vv = NULL;

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
    const PetscPointFn* strainKernel =
        (3 == spaceDim) ? pylith::fekernels::Elasticity3D::infinitesimalStrain_asVector :
        (2 == spaceDim) ? pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain_asVector :
        NULL;
    kernels[1] = ProjectKernels("cauchy_strain", strainKernel);

    const PetscPointFn* bulkDensity = pylith::fekernels::Poroelasticity::bulkDensity_asScalar;

    kernels[2] = ProjectKernels("bulk_density", bulkDensity);

    kernels[3] = ProjectKernels("water_content", _rheology->getKernelWaterContent(coordsys));

    assert(integrator);
    integrator->setKernelsDerivedField(kernels);

    PYLITH_METHOD_END;
} // _setKernelsDerivedField


// End of file
