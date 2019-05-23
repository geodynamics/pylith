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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesiveKin.hh" // implementation of object methods

#include "pylith/faults/KinSrc.hh" // USES KinSrc
#include "pylith/faults/AuxiliaryFactoryKinematic.hh" // USES AuxiliaryFactoryKinematic
#include "pylith/feassemble/IntegratorInterface.hh" // USES IntegratorInterface

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps::checkDiscretization()

#include "pylith/fekernels/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <cmath> // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring> // USES strlen()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error
#include <typeinfo> // USES typeid()

// #define DETAILED_EVENT_LOGGING

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorInterface::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorInterface::JacobianKernels JacobianKernels;

// ---------------------------------------------------------------------------------------------------------------------
namespace pylith {
    namespace faults {
        class _FaultCohesiveKin {
            // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /** Set kernels for RHS residual.
             *
             * @param[out] integrator Integrator for interface.
             * @param[in] fault Fault object for kinematic ruptures.
             * @param[in] solution Solution field.
             */
            static
            void setKernelsRHSResidual(pylith::feassemble::IntegratorInterface* integrator,
                                       const pylith::faults::FaultCohesiveKin& fault,
                                       const pylith::topology::Field& solution);

            /** Set kernels for RHS Jacobian.
             *
             * @param[out] integrator Integrator for interface.
             * @param[in] fault Fault object for kinematic ruptures.
             * @param[in] solution Solution field.
             */
            static
            void setKernelsRHSJacobian(pylith::feassemble::IntegratorInterface* integrator,
                                       const pylith::faults::FaultCohesiveKin& fault,
                                       const pylith::topology::Field& solution);

        };

        // _FaultCohesiveKin

    } // faults
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveKin::FaultCohesiveKin(void) :
    _auxiliaryFactory(new pylith::faults::AuxiliaryFactoryKinematic) {
    pylith::utils::PyreComponent::setName("faultcohesivekin");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveKin::~FaultCohesiveKin(void) { // destructor
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::FaultCohesiveKin::deallocate(void) {
    FaultCohesive::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
    _ruptures.clear(); // :TODO: Use shared pointers for earthquake ruptures
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set kinematic earthquake ruptures.
void
pylith::faults::FaultCohesiveKin::setEqRuptures(const char* const * names,
                                                const int numNames,
                                                KinSrc** ruptures,
                                                const int numRuptures) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setEqRuptures(names="<<names<<", numNames="<<numNames<<", ruptures="<<ruptures<<", numRuptures="<<numRuptures<<")");

    assert(numNames == numRuptures);

    // :TODO: Use shared pointers for earthquake ruptures
    _ruptures.clear();
    for (int i = 0; i < numRuptures; ++i) {
        if (!ruptures[i]) {
            std::ostringstream msg;
            msg << "Null earthquake rupture object for earthquake rupture '" << names[i] << "'.";
            throw std::runtime_error(msg.str());
        } // if
        _ruptures[std::string(names[i])] = ruptures[i];
    } // for

    PYLITH_METHOD_END;
} // setEqRuptures


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveKin::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    if (!solution.hasSubfield("displacement")) {
        std::ostringstream msg;
        msg << "Cannot find 'displacement' subfield in solution field for fault implementation in component '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    if (!solution.hasSubfield("lagrange_multiplier_fault")) {
        std::ostringstream msg;
        msg << "Cannot find 'lagrange_multiplier_fault' subfield in solution field for fault implementation in component '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::faults::FaultCohesiveKin::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.label()<<")");

    pylith::feassemble::IntegratorInterface* integrator = new pylith::feassemble::IntegratorInterface(this);assert(integrator);
    integrator->setInterfaceId(getInterfaceId());
    integrator->setSurfaceMarkerLabel(getSurfaceMarkerLabel());

    _FaultCohesiveKin::setKernelsRHSResidual(integrator, *this, solution);
    _FaultCohesiveKin::setKernelsRHSJacobian(integrator, *this, solution);
    // No state variables.
    // _FaultCohesiveKin::setKernelsDerivedFields(integrator, *this, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
pylith::feassemble::Constraint*
pylith::faults::FaultCohesiveKin::createConstraint(const pylith::topology::Field& solution) {
    return NULL;
} // createConstraint


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::faults::FaultCohesiveKin::createAuxiliaryField(const pylith::topology::Field& solution,
                                                       const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.label()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    assert(_normalizer);

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->label("FaultCohesiveKin auxiliary field");

    // Set default discretization of auxiliary subfields to match lagrange_multiplier_fault subfield in solution.
    assert(_auxiliaryFactory);
    const pylith::topology::FieldBase::Discretization& discretization = solution.subfieldInfo("lagrange_multiplier_fault").fe;
    _auxiliaryFactory->setSubfieldDiscretization("default", discretization.basisOrder, discretization.quadOrder, -1,
                                                 discretization.isBasisContinuous, discretization.feSpace);

    assert(_auxiliaryFactory);
    assert(_normalizer);
    _auxiliaryFactory->initialize(auxiliaryField, *_normalizer, domainMesh.dimension());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    _auxiliaryFactory->addSlip(); // 0

    auxiliaryField->subfieldsSetup();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->zeroLocal();

    // We don't populate the auxiliary field via a spatial database, because they will be set from the earthquake
    // rupture.

    // Initialize auxiliary fields for kinematic ruptures.
    assert(auxiliaryField);
    const srcs_type::const_iterator rupturesEnd = _ruptures.end();
    for (srcs_type::iterator r_iter = _ruptures.begin(); r_iter != rupturesEnd; ++r_iter) {
        KinSrc* src = r_iter->second;
        assert(src);
        src->initialize(*auxiliaryField, *_normalizer, solution.mesh().coordsys());
    } // for

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::faults::FaultCohesiveKin::createDerivedField(const pylith::topology::Field& solution,
                                                     const pylith::topology::Mesh& domainMesh) {
    return NULL;
} // createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::faults::FaultCohesiveKin::prestep(pylith::topology::Field* auxiliaryField,
                                          const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("prestep(auxiliaryField="<<auxiliaryField<<", t="<<t<<")");

    // Update slip subfield.

    // Compute slip field at current time step
    const srcs_type::const_iterator rupturesEnd = _ruptures.end();
    for (srcs_type::iterator r_iter = _ruptures.begin(); r_iter != rupturesEnd; ++r_iter) {
        KinSrc* src = r_iter->second;assert(src);
        src->slip(auxiliaryField, t);
    } // for

    PYLITH_METHOD_END;
} // prestep


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::faults::FaultCohesiveKin::_getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // _getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::faults::FaultCohesiveKin::_updateKernelConstants(const PylithReal dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelConstants(dt="<<dt<<")");

    if (6 != _kernelConstants.size()) { _kernelConstants.resize(6);}
    _kernelConstants[0] = _refDir1[0];
    _kernelConstants[1] = _refDir1[1];
    _kernelConstants[2] = _refDir1[2];
    _kernelConstants[3] = _refDir2[0];
    _kernelConstants[4] = _refDir2[1];
    _kernelConstants[5] = _refDir2[2];

    PYLITH_METHOD_END;
} // _updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS residual.
void
pylith::faults::_FaultCohesiveKin::setKernelsRHSResidual(pylith::feassemble::IntegratorInterface* integrator,
                                                         const pylith::faults::FaultCohesiveKin& fault,
                                                         const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    // PYLITH_COMPONENT_DEBUG("setKernelsRHSResidual(integrator="<<integrator<<", fault="<<typeid(fault)<<",
    // solution="<<solution.label()<<")");

    const char* nameDispVel = (solution.hasSubfield("velocity")) ? "velocity" : "displacement";
    const char* nameLagrangeMultiplier = "lagrange_multiplier_fault";

    std::vector<ResidualKernels> kernels(2);

    // Elasticity equation (displacement/velocity).
    const PetscBdPointFunc g0u = pylith::fekernels::FaultCohesiveKin::g0u;
    const PetscBdPointFunc g1u = NULL;

    // Fault slip constraint equation.
    const PetscBdPointFunc g0l = pylith::fekernels::FaultCohesiveKin::g0l;
    const PetscBdPointFunc g1l = NULL;

    kernels[0] = ResidualKernels(nameDispVel, g0u, g1u);
    kernels[1] = ResidualKernels(nameLagrangeMultiplier, g0l, g1l);
    assert(integrator);
    integrator->setKernelsRHSResidual(kernels);

    PYLITH_METHOD_END;
} // setKernelsRHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS Jacobian.
void
pylith::faults::_FaultCohesiveKin::setKernelsRHSJacobian(pylith::feassemble::IntegratorInterface* integrator,
                                                         const pylith::faults::FaultCohesiveKin& fault,
                                                         const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    // PYLITH_COMPONENT_DEBUG("setKernelsRHSJacobian(integrator="<<integrator<<", fault="<<typeid(fault)<<",
    // solution="<<solution.label()<<")");

    const char* nameDispVel = (solution.hasSubfield("velocity")) ? "velocity" : "displacement";
    const char* nameLagrangeMultiplier = "lagrange_multiplier_fault";

    std::vector<JacobianKernels> kernels(2);

    const PetscBdPointJac Jg0ul = pylith::fekernels::FaultCohesiveKin::Jg0ul;
    const PetscBdPointJac Jg1ul = NULL;
    const PetscBdPointJac Jg2ul = NULL;
    const PetscBdPointJac Jg3ul = NULL;

    const PetscBdPointJac Jg0lu = pylith::fekernels::FaultCohesiveKin::Jg0lu;
    const PetscBdPointJac Jg1lu = NULL;
    const PetscBdPointJac Jg2lu = NULL;
    const PetscBdPointJac Jg3lu = NULL;

    kernels[0] = JacobianKernels(nameDispVel, nameLagrangeMultiplier, Jg0ul, Jg1ul, Jg2ul, Jg3ul);
    kernels[1] = JacobianKernels(nameLagrangeMultiplier, nameDispVel, Jg0lu, Jg1lu, Jg2lu, Jg3lu);
    assert(integrator);
    integrator->setKernelsRHSJacobian(kernels);

    PYLITH_METHOD_END;
} // setKernelsRHSJacobian


// ----------------------------------------------------------------------
// End of file
