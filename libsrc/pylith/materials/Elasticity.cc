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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/materials/Elasticity.hh" // implementation of object methods

#include "pylith/materials/RheologyElasticity.hh" // HASA RheologyElasticity
#include "pylith/materials/AuxiliaryFactoryElasticity.hh" // USES AuxiliaryFactoryElasticity
#include "pylith/materials/DerivedFactoryElasticity.hh" // USES DerivedFactoryElasticity
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/feassemble/IntegratorInterface.hh" // USES IntegratorInterface
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldOps.hh" // USES FieldOps

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
typedef pylith::feassemble::Integrator::JacobianPart JacobianPart;
typedef pylith::feassemble::Integrator::ResidualPart ResidualPart;
typedef pylith::feassemble::IntegratorInterface::ResidualKernels InterfaceResidualKernels;
typedef pylith::feassemble::IntegratorInterface::JacobianKernels InterfaceJacobianKernels;
typedef pylith::feassemble::IntegratorInterface::FaceEnum InterfaceFace;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::Elasticity::Elasticity(void) :
    _useBodyForce(false),
    _rheology(NULL),
    _derivedFactory(new pylith::materials::DerivedFactoryElasticity) {
    pylith::utils::PyreComponent::setName("elasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::Elasticity::~Elasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::Elasticity::deallocate(void) {
    Material::deallocate();

    delete _derivedFactory;_derivedFactory = NULL;
    _rheology = NULL; // :TODO: Use shared pointer.
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
void
pylith::materials::Elasticity::useBodyForce(const bool value) {
    PYLITH_COMPONENT_DEBUG("useBodyForce(value="<<value<<")");

    _useBodyForce = value;
} // useBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
bool
pylith::materials::Elasticity::useBodyForce(void) const {
    return _useBodyForce;
} // useBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Set bulk rheology.
void
pylith::materials::Elasticity::setBulkRheology(pylith::materials::RheologyElasticity* const rheology) {
    _rheology = rheology;
} // setBulkRheology


// ---------------------------------------------------------------------------------------------------------------------
// Get bulk rheology.
pylith::materials::RheologyElasticity*
pylith::materials::Elasticity::getBulkRheology(void) const {
    return _rheology;
} // getBulkRheology


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::materials::Elasticity::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    // Verify solution contains expected fields.
    if (!solution.hasSubfield("displacement")) {
        throw std::runtime_error("Cannot find 'displacement' field in solution; required for 'Elasticity'.");
    } // if

    switch (_formulation) {
    case QUASISTATIC:
        break;
    case DYNAMIC:
    case DYNAMIC_IMEX:
        if (!solution.hasSubfield("velocity")) {
            throw std::runtime_error("Cannot find 'velocity' field in solution; required for 'dynamic' and "
                                     "'dynamic_imex' time stepping formulations of 'Elasticity'.");
        } // if
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::materials::Elasticity::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    pylith::feassemble::IntegratorDomain* integrator = new pylith::feassemble::IntegratorDomain(this);assert(integrator);
    integrator->setLabelName(pylith::topology::Mesh::getCellsLabelName());
    integrator->setLabelValue(getMaterialId());

    _setKernelsResidual(integrator, solution);
    _setKernelsJacobian(integrator, solution);
    _setKernelsUpdateStateVars(integrator, solution);
    _setKernelsDerivedField(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::materials::Elasticity::createAuxiliaryField(const pylith::topology::Field& solution,
                                                    const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getLabel()<<", domainMesh="<<typeid(domainMesh).name()<<")");

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->setLabel("Elasticity auxiliary field");

    assert(_rheology);
    pylith::materials::AuxiliaryFactoryElasticity* auxiliaryFactory = _rheology->getAuxiliaryFactory();assert(auxiliaryFactory);

    assert(_normalizer);
    auxiliaryFactory->initialize(auxiliaryField, *_normalizer, domainMesh.getDimension());

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
pylith::materials::Elasticity::createDerivedField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.getLabel()<<", domainMesh="<<typeid(domainMesh).name()<<")");

    assert(_derivedFactory);
    if (_derivedFactory->getNumSubfields() == 1) {
        PYLITH_METHOD_RETURN(NULL);
    } // if

    pylith::topology::Field* derivedField = new pylith::topology::Field(domainMesh);assert(derivedField);
    derivedField->setLabel("Elasticity derived field");

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


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::materials::Elasticity::_getAuxiliaryFactory(void) {
    assert(_rheology);
    return _rheology->getAuxiliaryFactory();
} // _getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::Elasticity::_updateKernelConstants(const PylithReal dt) {
    assert(_rheology);
    _rheology->updateKernelConstants(&_kernelConstants, dt);
} // _updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Get derived factory associated with physics.
pylith::topology::FieldFactory*
pylith::materials::Elasticity::_getDerivedFactory(void) {
    return _derivedFactory;
} // _getDerivedFactory


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::materials::Elasticity::_setKernelsResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                   const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsResidual(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();

    const int bitBodyForce = _useBodyForce ? 0x1 : 0x0;
    const int bitGravity = _gravityField ? 0x2 : 0x0;
    const int bitUse = bitBodyForce | bitGravity;

    PetscPointFunc r0 = NULL;
    PetscBdPointFunc l0_neg_bodyforce = NULL;
    PetscBdPointFunc l0_pos_bodyforce = NULL;
    switch (bitUse) {
    case 0x1:
        r0 = pylith::fekernels::Elasticity::g0v_bodyforce;
        l0_neg_bodyforce = pylith::fekernels::Elasticity::f0l_neg_bodyforce;
        l0_pos_bodyforce = pylith::fekernels::Elasticity::f0l_pos_bodyforce;
        break;
    case 0x2:
        r0 = pylith::fekernels::Elasticity::g0v_grav;
        l0_neg_bodyforce = pylith::fekernels::Elasticity::f0l_neg_grav;
        l0_pos_bodyforce = pylith::fekernels::Elasticity::f0l_pos_grav;
        break;
    case 0x3:
        r0 = pylith::fekernels::Elasticity::g0v_gravbodyforce;
        l0_neg_bodyforce = pylith::fekernels::Elasticity::f0l_neg_gravbodyforce;
        l0_pos_bodyforce = pylith::fekernels::Elasticity::f0l_pos_gravbodyforce;
        break;
    case 0x0:
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown case (bitUse=" << bitUse << ") for residual kernels.");
    } // switch
    const PetscPointFunc r1 = _rheology->getKernelResidualStress(coordsys);

    std::vector<ResidualKernels> kernels;
    std::vector<InterfaceResidualKernels> kernelsInterface;
    switch (_formulation) {
    case QUASISTATIC: {
        const PetscPointFunc f0u = r0;
        const PetscPointFunc f1u = r1;

        kernels.resize(1);
        kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::RESIDUAL_LHS, f0u, f1u);
        break;
    } // QUASISTATIC
    case DYNAMIC_IMEX: {
        // Displacement
        const PetscPointFunc f0u = pylith::fekernels::DispVel::f0u;
        const PetscPointFunc f1u = NULL;
        const PetscPointFunc g0u = pylith::fekernels::DispVel::g0u;
        const PetscPointFunc g1u = NULL;

        // Velocity
        const PetscPointFunc f0v = pylith::fekernels::DispVel::f0v;
        const PetscPointFunc f1v = NULL;
        const PetscPointFunc g0v = r0;
        const PetscPointFunc g1v = r1;

        kernels.resize(4);
        kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::RESIDUAL_LHS, f0u, f1u);
        kernels[1] = ResidualKernels("displacement", pylith::feassemble::Integrator::RESIDUAL_RHS, g0u, g1u);
        kernels[2] = ResidualKernels("velocity", pylith::feassemble::Integrator::RESIDUAL_LHS, f0v, f1v);
        kernels[3] = ResidualKernels("velocity", pylith::feassemble::Integrator::RESIDUAL_RHS, g0v, g1v);

        // Interface kernels
        PetscBdPointFunc f0l_neg = pylith::fekernels::Elasticity::f0l_neg;
        PetscBdPointFunc f1l_neg = NULL;
        PetscBdPointFunc f0l_pos = pylith::fekernels::Elasticity::f0l_pos;
        PetscBdPointFunc f1l_pos = NULL;

        const PetscBdPointFunc f0l_stress_neg = _rheology->getInterfaceKernelResidualF0Neg(coordsys);
        const PetscBdPointFunc f1l_stress_neg = _rheology->getInterfaceKernelResidualF1Neg(coordsys);
        const PetscBdPointFunc f0l_stress_pos = _rheology->getInterfaceKernelResidualF0Pos(coordsys);
        const PetscBdPointFunc f1l_stress_pos = _rheology->getInterfaceKernelResidualF1Pos(coordsys);

        kernelsInterface.resize(4);
        const ResidualPart residualTerm = pylith::feassemble::Integrator::RESIDUAL_LHS;
        const InterfaceFace faceNeg = pylith::feassemble::IntegratorInterface::NEGATIVE_FACE;
        const InterfaceFace facePos = pylith::feassemble::IntegratorInterface::POSITIVE_FACE;
        kernelsInterface[0] = InterfaceResidualKernels("lagrange_multiplier_fault", residualTerm, faceNeg, f0l_neg, f1l_neg);
        kernelsInterface[1] = InterfaceResidualKernels("lagrange_multiplier_fault", residualTerm, facePos, f0l_pos, f1l_pos);
        kernelsInterface[2] = InterfaceResidualKernels("lagrange_multiplier_fault", residualTerm, faceNeg, f0l_stress_neg, f1l_stress_neg);
        kernelsInterface[3] = InterfaceResidualKernels("lagrange_multiplier_fault", residualTerm, facePos, f0l_stress_pos, f1l_stress_pos);
        break;
    } // DYNAMIC_IMEX
    case DYNAMIC: {
        // Displacement
        const PetscPointFunc f0u = pylith::fekernels::DispVel::f0u;
        const PetscPointFunc f1u = NULL;
        const PetscPointFunc g0u = pylith::fekernels::DispVel::g0u;
        const PetscPointFunc g1u = NULL;

        // Velocity
        const PetscPointFunc f0v = pylith::fekernels::DispVel::f0v;
        const PetscPointFunc f1v = NULL;
        const PetscPointFunc g0v = r0;
        const PetscPointFunc g1v = r1;

        kernels.resize(4);
        kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::RESIDUAL_LHS, f0u, f1u);
        kernels[1] = ResidualKernels("displacement", pylith::feassemble::Integrator::RESIDUAL_RHS, g0u, g1u);
        kernels[2] = ResidualKernels("velocity", pylith::feassemble::Integrator::RESIDUAL_LHS, f0v, f1v);
        kernels[3] = ResidualKernels("velocity", pylith::feassemble::Integrator::RESIDUAL_RHS, g0v, g1v);
        break;
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);
    integrator->setKernelsResidual(kernelsInterface, solution);

    PYLITH_METHOD_END;
} // _setKernelsResidual


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for Jacobian.
void
pylith::materials::Elasticity::_setKernelsJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                   const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsJacobian(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");
    assert(integrator);

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();

    std::vector<JacobianKernels> kernels;
    std::vector<InterfaceJacobianKernels> kernelsInterface;

    switch (_formulation) {
    case QUASISTATIC: {
        const PetscPointJac Jf0uu = NULL;
        const PetscPointJac Jf1uu = NULL;
        const PetscPointJac Jf2uu = NULL;
        const PetscPointJac Jf3uu = _rheology->getKernelJacobianElasticConstants(coordsys);

        integrator->setLHSJacobianTriggers(_rheology->getLHSJacobianTriggers());

        kernels.resize(1);
        const JacobianPart jacobianPart = pylith::feassemble::Integrator::JACOBIAN_LHS;
        kernels[0] = JacobianKernels("displacement", "displacement", jacobianPart, Jf0uu, Jf1uu, Jf2uu, Jf3uu);
        break;
    } // QUASISTATIC
    case DYNAMIC_IMEX: {
        kernelsInterface.resize(2);

        const PetscBdPointJac Jf0lu_neg = NULL;
        const PetscBdPointJac Jf1lu_neg = _rheology->getInterfaceKernelJacobianF1Neg(coordsys);
        const PetscBdPointJac Jf2lu_neg = NULL;
        const PetscBdPointJac Jf3lu_neg = _rheology->getInterfaceKernelJacobianF3Neg(coordsys);

        const PetscBdPointJac Jf0lu_pos = NULL;
        const PetscBdPointJac Jf1lu_pos = _rheology->getInterfaceKernelJacobianF1Pos(coordsys);
        const PetscBdPointJac Jf2lu_pos = NULL;
        const PetscBdPointJac Jf3lu_pos = _rheology->getInterfaceKernelJacobianF3Pos(coordsys);

        kernelsInterface.resize(2);
        const JacobianPart jacobianPart = pylith::feassemble::Integrator::JACOBIAN_LHS;
        const InterfaceFace faceNeg = pylith::feassemble::IntegratorInterface::NEGATIVE_FACE;
        const InterfaceFace facePos = pylith::feassemble::IntegratorInterface::POSITIVE_FACE;
        kernelsInterface[0] = InterfaceJacobianKernels("lagrange_multiplier_fault", "displacement", jacobianPart, faceNeg,
                                                       Jf0lu_neg, Jf1lu_neg, Jf2lu_neg, Jf3lu_neg);
        kernelsInterface[1] = InterfaceJacobianKernels("lagrange_multiplier_fault", "displacement", jacobianPart, facePos,
                                                       Jf0lu_pos, Jf1lu_pos, Jf2lu_pos, Jf3lu_pos);
    } // DYNAMIC_IMEX (continue with DYNAMIC)
    case DYNAMIC: {
        const PetscPointJac Jf0uu = pylith::fekernels::DispVel::Jf0uu_stshift;
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

        const PetscPointJac Jf0vv = pylith::fekernels::Elasticity::Jf0vv;
        const PetscPointJac Jf1vv = NULL;
        const PetscPointJac Jf2vv = NULL;
        const PetscPointJac Jf3vv = NULL;

        integrator->setLHSJacobianTriggers(pylith::feassemble::Integrator::NEW_JACOBIAN_TIME_STEP_CHANGE);

        kernels.resize(4);
        const JacobianPart jacobianPart = pylith::feassemble::Integrator::JACOBIAN_LHS_LUMPED_INV;
        kernels[0] = JacobianKernels("displacement", "displacement", jacobianPart, Jf0uu, Jf1uu, Jf2uu, Jf3uu);
        kernels[1] = JacobianKernels("displacement", "velocity", jacobianPart, Jf0uv, Jf1uv, Jf2uv, Jf3uv);
        kernels[2] = JacobianKernels("velocity", "displacement", jacobianPart, Jf0vu, Jf1vu, Jf2vu, Jf3vu);
        kernels[3] = JacobianKernels("velocity", "velocity", jacobianPart, Jf0vv, Jf1vv, Jf2vv, Jf3vv);
        break;
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    integrator->setKernelsJacobian(kernels, solution);
    integrator->setKernelsJacobian(kernelsInterface, solution);

    PYLITH_METHOD_END;
} // _setKernelsJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for computing updated state variables in auxiliary field.
void
pylith::materials::Elasticity::_setKernelsUpdateStateVars(pylith::feassemble::IntegratorDomain* integrator,
                                                          const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsUpdateStateVars(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    std::vector<ProjectKernels> kernels;
    _rheology->addKernelsUpdateStateVars(&kernels, coordsys);

    integrator->setKernelsUpdateStateVars(kernels);

    PYLITH_METHOD_END;
} // _setKernelsUpdateStateVars


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for computing derived field.
void
pylith::materials::Elasticity::_setKernelsDerivedField(pylith::feassemble::IntegratorDomain* integrator,
                                                       const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsDerivedField(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    std::vector<ProjectKernels> kernels(2);
    kernels[0] = ProjectKernels("cauchy_stress", _rheology->getKernelDerivedCauchyStress(coordsys));

    const int spaceDim = coordsys->getSpaceDim();
    const PetscPointFunc strainKernel =
        (3 == spaceDim) ? pylith::fekernels::Elasticity3D::cauchyStrain :
        (2 == spaceDim) ? pylith::fekernels::ElasticityPlaneStrain::cauchyStrain :
        NULL;
    kernels[1] = ProjectKernels("cauchy_strain", strainKernel);

    assert(integrator);
    integrator->setKernelsDerivedField(kernels);

    PYLITH_METHOD_END;
} // _setKernelsDerivedField


// End of file
