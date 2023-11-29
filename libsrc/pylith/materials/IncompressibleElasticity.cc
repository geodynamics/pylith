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

#include "pylith/materials/IncompressibleElasticity.hh" // implementation of object methods

#include "pylith/materials/RheologyIncompressibleElasticity.hh" // HASA RheologyIncompressibleElasticity
#include "pylith/materials/AuxiliaryFactoryElasticity.hh" // USES AuxiliaryFactoryElasticity
#include "pylith/materials/DerivedFactoryElasticity.hh" // USES DerivedFactoryElasticity
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/fekernels/IncompressibleElasticity.hh" // USES IncompressibleElasticity kernels
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/DispVel.hh" // USES DispVel kernels

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;
typedef pylith::feassemble::Integrator::EquationPart EquationPart;

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IncompressibleElasticity::IncompressibleElasticity(void) :
    _useBodyForce(false),
    _rheology(NULL),
    _derivedFactory(new pylith::materials::DerivedFactoryElasticity) {
    pylith::utils::PyreComponent::setName("incompressibleelasticity");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IncompressibleElasticity::~IncompressibleElasticity(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IncompressibleElasticity::deallocate(void) {
    Material::deallocate();

    delete _derivedFactory;_derivedFactory = NULL;
    _rheology = NULL; // :TODO: Use shared pointer.
} // deallocate


// ------------------------------------------------------------------------------------------------
// Include body force?
void
pylith::materials::IncompressibleElasticity::useBodyForce(const bool value) {
    PYLITH_COMPONENT_DEBUG("useBodyForce(value="<<value<<")");

    _useBodyForce = value;
} // useBodyForce


// ------------------------------------------------------------------------------------------------
// Include body force?
bool
pylith::materials::IncompressibleElasticity::useBodyForce(void) const {
    return _useBodyForce;
} // useBodyForce


// ------------------------------------------------------------------------------------------------
// Set bulk rheology.
void
pylith::materials::IncompressibleElasticity::setBulkRheology(pylith::materials::RheologyIncompressibleElasticity* const rheology) {
    PYLITH_COMPONENT_DEBUG("setBulkRheology(rheology="<<rheology<<")");

    _rheology = rheology;
} // setBulkRheology


// ------------------------------------------------------------------------------------------------
// Get bulk rheology.
pylith::materials::RheologyIncompressibleElasticity*
pylith::materials::IncompressibleElasticity::getBulkRheology(void) const {
    return _rheology;
} // getBulkRheology


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::materials::IncompressibleElasticity::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    // Verify solution contains required fields.
    std::string reason = "material 'IncompressibleElasticity'.";
    size_t numRequired = 0;
    const size_t maxRequired = 2;
    pylith::string_vector requiredFields(maxRequired);
    requiredFields[numRequired++] = "displacement";
    requiredFields[numRequired++] = "pressure";
    requiredFields.resize(numRequired);

    pylith::topology::FieldOps::checkSubfieldsExist(requiredFields, reason, solution);

    PYLITH_METHOD_END;
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::materials::IncompressibleElasticity::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    pylith::feassemble::IntegratorDomain* integrator = new pylith::feassemble::IntegratorDomain(this);assert(integrator);
    integrator->setLabelName(getLabelName());
    integrator->setLabelValue(getLabelValue());

    _setKernelsResidual(integrator, solution);
    _setKernelsJacobian(integrator, solution);
    // No state variables.
    _setKernelsDerivedField(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::materials::IncompressibleElasticity::createAuxiliaryField(const pylith::topology::Field& solution,
                                                                  const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->setLabel("IncompressibleElasticity auxiliary field");

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


// ------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::materials::IncompressibleElasticity::createDerivedField(const pylith::topology::Field& solution,
                                                                const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

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


// ------------------------------------------------------------------------------------------------
// Get default PETSc solver options appropriate for material.
pylith::utils::PetscOptions*
pylith::materials::IncompressibleElasticity::getSolverDefaults(const bool isParallel,
                                                               const bool hasFault) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getSolverDefaults(isParallel="<<isParallel<<", hasFault="<<hasFault<<")");

    pylith::utils::PetscOptions* options = new pylith::utils::PetscOptions();assert(options);

    switch (_formulation) {
    case pylith::problems::Physics::QUASISTATIC:
        options->add("-ts_type", "beuler");

        if (!hasFault) {
            options->add("-pc_type", "fieldsplit");
            options->add("-pc_fieldsplit_type", "schur");
            options->add("-pc_fieldsplit_schur_factorization_type", "full");
            options->add("-pc_fieldsplit_schur_precondition", "full");
            if (!isParallel) {
                options->add("-fieldsplit_displacement_pc_type", "lu");
                options->add("-fieldsplit_pressure_pc_type", "lu");
            } else {
#if 1
                options->add("-fieldsplit_displacement_pc_type", "ml");
#else
                options->add("-fieldsplit_displacement_pc_type", "gamg");
                options->add("-fieldsplit_displacement_mg_levels_pc_type", "sor");
                options->add("-fieldsplit_displacement_mg_levels_ksp_type", "richardson");
#endif
                options->add("-fieldsplit_pressure_pc_type", "bjacobi");
            } // if/else
        } // if/else
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


// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::materials::IncompressibleElasticity::_getAuxiliaryFactory(void) {
    assert(_rheology);
    return _rheology->getAuxiliaryFactory();
} // _getAuxiliaryFactory


// ------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::IncompressibleElasticity::_updateKernelConstants(const PylithReal dt) {
    assert(_rheology);
    _rheology->updateKernelConstants(&_kernelConstants, dt);
} // _updateKernelConstants


// ------------------------------------------------------------------------------------------------
// Get derived factory associated with physics.
pylith::topology::FieldFactory*
pylith::materials::IncompressibleElasticity::_getDerivedFactory(void) {
    return _derivedFactory;
} // _getDerivedFactory


// ------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::materials::IncompressibleElasticity::_setKernelsResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                                 const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsResidual(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();

    const int bitBodyForce = _useBodyForce ? 0x1 : 0x0;
    const int bitGravity = _gravityField ? 0x2 : 0x0;
    const int bitUse = bitBodyForce | bitGravity;

    PetscPointFunc r0 = NULL;
    switch (bitUse) {
    case 0x1:
        r0 = pylith::fekernels::Elasticity::g0v_bodyforce;
        break;
    case 0x2:
        r0 = pylith::fekernels::Elasticity::g0v_grav;
        break;
    case 0x3:
        r0 = pylith::fekernels::Elasticity::g0v_gravbodyforce;
        break;
    case 0x0:
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown case (bitUse=" << bitUse << ") for residual kernels.");
    } // switch

    // Displacement
    const PetscPointFunc f0u = r0;
    const PetscPointFunc f1u = _rheology->getKernelf1u(coordsys);

    // Pressure
    const PetscPointFunc f0p = _rheology->getKernelf0p(coordsys);
    const PetscPointFunc f1p = NULL;

    std::vector<ResidualKernels> kernels(2);
    kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::LHS, f0u, f1u);
    kernels[1] = ResidualKernels("pressure", pylith::feassemble::Integrator::LHS, f0p, f1p);

    // Add any MMS body force kernels.
    kernels.insert(kernels.end(), _mmsBodyForceKernels.begin(), _mmsBodyForceKernels.end());

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsResidual


// ------------------------------------------------------------------------------------------------
// Set kernels for Jacobian.
void
pylith::materials::IncompressibleElasticity::_setKernelsJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                                 const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsJacobian(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();

    std::vector<JacobianKernels> kernels(4);

    const PetscPointJac Jf0uu = NULL;
    const PetscPointJac Jf1uu = NULL;
    const PetscPointJac Jf2uu = NULL;
    const PetscPointJac Jf3uu = _rheology->getKernelJf3uu(coordsys);

    const PetscPointJac Jf0up = NULL;
    const PetscPointJac Jf1up = NULL;
    const PetscPointJac Jf2up = pylith::fekernels::IncompressibleElasticity::Jf2up;
    const PetscPointJac Jf3up = NULL;

    const PetscPointJac Jf0pu = NULL;
    const PetscPointJac Jf1pu = pylith::fekernels::IncompressibleElasticity::Jf1pu;
    const PetscPointJac Jf2pu = NULL;
    const PetscPointJac Jf3pu = NULL;

    const PetscPointJac Jf0pp = _rheology->getKernelJf0pp(coordsys);
    const PetscPointJac Jf1pp = NULL;
    const PetscPointJac Jf2pp = NULL;
    const PetscPointJac Jf3pp = NULL;

    const EquationPart equationPart = pylith::feassemble::Integrator::LHS;
    kernels[0] = JacobianKernels("displacement", "displacement", equationPart, Jf0uu, Jf1uu, Jf2uu, Jf3uu);
    kernels[1] = JacobianKernels("displacement", "pressure", equationPart, Jf0up, Jf1up, Jf2up, Jf3up);
    kernels[2] = JacobianKernels("pressure", "displacement", equationPart, Jf0pu, Jf1pu, Jf2pu, Jf3pu);
    kernels[3] = JacobianKernels("pressure", "pressure", equationPart, Jf0pp, Jf1pp, Jf2pp, Jf3pp);

    assert(integrator);
    integrator->setKernelsJacobian(kernels, solution);

    PYLITH_METHOD_END;
} // setKernelsJacobian


// ------------------------------------------------------------------------------------------------
// Set kernels for computing updated state variables in auxiliary field.
void
pylith::materials::IncompressibleElasticity::_setKernelsUpdateStateVars(pylith::feassemble::IntegratorDomain* integrator,
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


// ------------------------------------------------------------------------------------------------
// Set kernels for computing derived field.
void
pylith::materials::IncompressibleElasticity::_setKernelsDerivedField(pylith::feassemble::IntegratorDomain* integrator,
                                                                     const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsDerivedField(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    std::vector<ProjectKernels> kernels(2);
    kernels[0] = ProjectKernels("cauchy_stress", _rheology->getKernelCauchyStressVector(coordsys));

    const int spaceDim = coordsys->getSpaceDim();
    const PetscPointFunc strainKernel =
        (3 == spaceDim) ? pylith::fekernels::Elasticity3D::infinitesimalStrain_asVector :
        (2 == spaceDim) ? pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain_asVector :
        NULL;
    kernels[1] = ProjectKernels("cauchy_strain", strainKernel);

    assert(integrator);
    integrator->setKernelsDerivedField(kernels);

    PYLITH_METHOD_END;
} // _setKernelsDerivedField


// End of file
