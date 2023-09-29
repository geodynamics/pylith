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

#include "FaultCohesiveImpulses.hh" // implementation of object methods

#include "pylith/faults/AuxiliaryFactoryKinematic.hh" // USES AuxiliaryFactoryKinematic
#include "pylith/feassemble/IntegratorInterface.hh" // USES IntegratorInterface
#include "pylith/feassemble/InterfacePatches.hh" // USES InterfacePatches
#include "pylith/feassemble/ConstraintSimple.hh" // USES ConstraintSimple

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps::checkDiscretization()
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/fekernels/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensionalizer
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <cmath> // USES pow(), sqrt()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error
#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorInterface::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorInterface::JacobianKernels JacobianKernels;

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace faults {
        class _FaultCohesiveImpulses {
            // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////
public:

            static const char* pyreComponent;

            /** Find points with impulse amplitude greater than threshold.
             *
             * @param[out] impulsePoints Array of points on which to apply impulses.
             * @param[in] auxiliaryField Auxiliary field with impulse amplitude.
             * @param[in] threshold Threshold for impulses.
             */
            static
            void findImpulsePoints(int_array* impulsePoints,
                                   const pylith::topology::Field& auxiliaryField,
                                   const double threshold);

        };
        const char* _FaultCohesiveImpulses::pyreComponent = "faultcohesiveimpulses";

        // _FaultCohesiveKin

    } // faults
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveImpulses::FaultCohesiveImpulses(void) :
    _auxiliaryFactory(new pylith::faults::AuxiliaryFactoryKinematic),
    _threshold(1.0e-6) {
    pylith::utils::PyreComponent::setName(_FaultCohesiveImpulses::pyreComponent);
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveImpulses::~FaultCohesiveImpulses(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::FaultCohesiveImpulses::deallocate(void) {
    FaultCohesive::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set indices of fault degrees of freedom associated with
void
pylith::faults::FaultCohesiveImpulses::setImpulseDOF(const int* flags,
                                                     const size_t size) {
    PYLITH_COMPONENT_DEBUG("setImpulseDOF(flags="<<flags<<", size="<<size<<")");

    if (((size == 0) && flags) || ((size > 0) && !flags)) {
        PYLITH_COMPONENT_LOGICERROR("Inconsistent array size ("<<size<<") and values ("<<flags<<").");
    } // if

    _impulseDOF.resize(size);
    for (size_t i = 0; i < size; ++i) {
        _impulseDOF[i] = flags[i];
    } // for
} // setImpulseDOF


// ------------------------------------------------------------------------------------------------
// Set threshold for nonzero impulse amplitude.
void
pylith::faults::FaultCohesiveImpulses::setThreshold(const double value) {
    PYLITH_COMPONENT_DEBUG("setThreshold(value="<<value<<")");

    if (value < 0) {
        std::ostringstream msg;
        msg << "Threshold ("<< value << ") for impulse amplitude must be nonnegative.";
        throw std::out_of_range(msg.str());
    } // if
    assert(_normalizer);
    _threshold = value / _normalizer->getLengthScale();
} // setThreshold


// ------------------------------------------------------------------------------------------------
// Return the number of impulses on this process.
size_t
pylith::faults::FaultCohesiveImpulses::getNumImpulsesLocal(void) {
    return _impulsePoints.size()*_impulseDOF.size();
} // getNumImpulses


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveImpulses::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    // Verify solution contains required fields.
    std::string reason = "interface 'FaultCohesiveImpulses'.";
    size_t numRequired = 0;
    const size_t maxRequired = 2;
    pylith::string_vector requiredFields(maxRequired);
    requiredFields[numRequired++] = "displacement";
    requiredFields[numRequired++] = "lagrange_multiplier_fault";
    requiredFields.resize(numRequired);

    pylith::topology::FieldOps::checkSubfieldsExist(requiredFields, reason, solution);

    PYLITH_METHOD_END;
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::faults::FaultCohesiveImpulses::createAuxiliaryField(const pylith::topology::Field& solution,
                                                            const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    assert(_normalizer);

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->setLabel("FaultCohesiveImpulses auxiliary field");

    // Set default discretization of auxiliary subfields to match lagrange_multiplier_fault subfield in solution.
    const pylith::topology::FieldBase::Discretization& discretization = solution.getSubfieldInfo("lagrange_multiplier_fault").fe;
    const PylithInt dimension = -1;
    const bool isFaultOnly = false;
    assert(_auxiliaryFactory);
    _auxiliaryFactory->setSubfieldDiscretization("default", discretization.basisOrder, discretization.quadOrder, dimension,
                                                 isFaultOnly, discretization.cellBasis, discretization.feSpace,
                                                 discretization.isBasisContinuous);

    assert(_normalizer);
    _auxiliaryFactory->initialize(auxiliaryField, *_normalizer, solution.getSpaceDim());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    switch (_formulation) {
    case QUASISTATIC:
        _auxiliaryFactory->addSlip(); // 0
        _auxiliaryFactory->setSubfieldQuery("slip"); // Not set in factory (assumes FaultCohesiveKin)
        break;
    case DYNAMIC_IMEX:
    case DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Green's functions for dynamic simulations is not yet supported.");
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();
    auxiliaryField->createGlobalVector();

    assert(_auxiliaryFactory);
    _auxiliaryFactory->setValuesFromDB();
    _FaultCohesiveImpulses::findImpulsePoints(&_impulsePoints, *auxiliaryField, _threshold);

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG("Displaying auxiliary field");
        auxiliaryField->view("FaultCohesiveImpulses auxiliary field");
    } // if

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::faults::FaultCohesiveImpulses::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                                            const double impulseReal) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateAuxiliaryField(auxiliaryField="<<auxiliaryField<<", impulseReal="<<impulseReal<<")");

    assert(auxiliaryField);
    assert(_normalizer);

    switch (_formulation) {
    case QUASISTATIC:
        this->_updateSlip(auxiliaryField, long(impulseReal));
        break;
    case DYNAMIC_IMEX:
    case DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Green's functions for dynamic simulations is not yet supported.");
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    PYLITH_METHOD_END;
} // updateAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::faults::FaultCohesiveImpulses::_getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // _getAuxiliaryFactory


// ------------------------------------------------------------------------------------------------
// Update slip subfield in auxiliary field at beginning of time step.
void
pylith::faults::FaultCohesiveImpulses::_updateSlip(pylith::topology::Field* auxiliaryField,
                                                   const long impulseStep) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_updateSlip(auxiliaryField="<<auxiliaryField<<", impulseStep="<<impulseStep<<")");

    assert(auxiliaryField);
    assert(_normalizer);

    const size_t numComponents = _impulseDOF.size();
    const size_t iComponent = impulseStep % numComponents;
    const size_t iImpulse = impulseStep / numComponents;

    PetscErrorCode err;
    err = VecSet(auxiliaryField->getLocalVector(), 0.0);PYLITH_CHECK_ERROR(err);

    if ((impulseStep >= 0) && (_impulsePoints.size() > 0)) {
        pylith::topology::VecVisitorMesh auxiliaryVisitor(*auxiliaryField, "slip");
        PylithScalar* auxiliaryArray = auxiliaryVisitor.localArray();
        const PetscInt pImpulse = _impulsePoints[iImpulse];
        const PetscInt dof = _impulseDOF[iComponent];
        const PetscInt slipDof = auxiliaryVisitor.sectionDof(pImpulse);assert(iComponent < size_t(slipDof));
        const PetscInt slipOff = auxiliaryVisitor.sectionOffset(pImpulse);
        auxiliaryArray[slipOff+dof] = 1.0 / _normalizer->getLengthScale();
    } // if

    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state()) {
        auxiliaryField->view("Fault auxiliary field after setting impulse.");
    } // if

    err = VecSet(auxiliaryField->getGlobalVector(), 0.0);PYLITH_CHECK_ERROR(err);
    auxiliaryField->scatterLocalToVector(auxiliaryField->getGlobalVector(), ADD_VALUES);
    auxiliaryField->scatterVectorToLocal(auxiliaryField->getGlobalVector(), INSERT_VALUES);

    PYLITH_METHOD_END;
} // _updateSlip


// ------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::faults::FaultCohesiveImpulses::_setKernelsResidual(pylith::feassemble::IntegratorInterface* integrator,
                                                           const pylith::topology::Field& solution,
                                                           const std::vector<pylith::materials::Material*>& materials) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsResidual(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");
    typedef pylith::feassemble::IntegratorInterface integrator_t;

    std::vector<ResidualKernels> kernels;
    switch (_formulation) {
    case pylith::problems::Physics::QUASISTATIC: {
        // Elasticity equation (displacement) for negative side of the fault.
        const PetscBdPointFunc f0u_neg = pylith::fekernels::FaultCohesiveKin::f0u_neg;
        const PetscBdPointFunc f1u_neg = NULL;

        // Elasticity equation (displacement) for positive side of the fault.
        const PetscBdPointFunc f0u_pos = pylith::fekernels::FaultCohesiveKin::f0u_pos;
        const PetscBdPointFunc f1u_pos = NULL;

        // Fault slip constraint equation.
        const PetscBdPointFunc f0l = pylith::fekernels::FaultCohesiveKin::f0l_slip;
        const PetscBdPointFunc f1l = NULL;

        kernels.resize(3);
        kernels[0] = ResidualKernels("displacement", integrator_t::LHS, integrator_t::NEGATIVE_FACE,
                                     f0u_neg, f1u_neg);
        kernels[1] = ResidualKernels("displacement", integrator_t::LHS, integrator_t::POSITIVE_FACE,
                                     f0u_pos, f1u_pos);
        kernels[2] = ResidualKernels("lagrange_multiplier_fault", integrator_t::LHS, integrator_t::FAULT_FACE,
                                     f0l, f1l);

        break;
    } // QUASISTATIC
    case pylith::problems::Physics::DYNAMIC_IMEX:
    case pylith::problems::Physics::DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Green's functions for dynamic simulations is not yet supported.");
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation '"<<_formulation<<"'.");
    } // switch

    assert(integrator);
    integrator->setKernels(kernels, solution, materials);

    PYLITH_METHOD_END;
} // _setKernelsResidual


// ------------------------------------------------------------------------------------------------
// Set kernels for Jacobian.
void
pylith::faults::FaultCohesiveImpulses::_setKernelsJacobian(pylith::feassemble::IntegratorInterface* integrator,
                                                           const pylith::topology::Field& solution,
                                                           const std::vector<pylith::materials::Material*>& materials) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsJacobian(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");
    typedef pylith::feassemble::IntegratorInterface integrator_t;

    std::vector<JacobianKernels> kernels;
    switch (_formulation) {
    case QUASISTATIC: {
        const PetscBdPointJac Jf0ul_neg = pylith::fekernels::FaultCohesiveKin::Jf0ul_neg;
        const PetscBdPointJac Jf1ul_neg = NULL;
        const PetscBdPointJac Jf2ul_neg = NULL;
        const PetscBdPointJac Jf3ul_neg = NULL;

        const PetscBdPointJac Jf0ul_pos = pylith::fekernels::FaultCohesiveKin::Jf0ul_pos;
        const PetscBdPointJac Jf1ul_pos = NULL;
        const PetscBdPointJac Jf2ul_pos = NULL;
        const PetscBdPointJac Jf3ul_pos = NULL;

        const PetscBdPointJac Jf0lu = pylith::fekernels::FaultCohesiveKin::Jf0lu;
        const PetscBdPointJac Jf1lu = NULL;
        const PetscBdPointJac Jf2lu = NULL;
        const PetscBdPointJac Jf3lu = NULL;

        kernels.resize(3);
        const char* nameDisplacement = "displacement";
        const char* nameLagrangeMultiplier = "lagrange_multiplier_fault";
        kernels[0] = JacobianKernels(nameDisplacement, nameLagrangeMultiplier, integrator_t::LHS,
                                     integrator_t::NEGATIVE_FACE, Jf0ul_neg, Jf1ul_neg, Jf2ul_neg, Jf3ul_neg);
        kernels[1] = JacobianKernels(nameDisplacement, nameLagrangeMultiplier, integrator_t::LHS,
                                     integrator_t::POSITIVE_FACE, Jf0ul_pos, Jf1ul_pos, Jf2ul_pos, Jf3ul_pos);
        kernels[2] = JacobianKernels(nameLagrangeMultiplier, nameDisplacement, integrator_t::LHS,
                                     integrator_t::FAULT_FACE, Jf0lu, Jf1lu, Jf2lu, Jf3lu);
        break;
    } // QUASISTATIC
    case pylith::problems::Physics::DYNAMIC_IMEX:
    case pylith::problems::Physics::DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Green's functions for dynamic simulations is not yet supported.");
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation '"<<_formulation<<"'.");
    } // switch

    assert(integrator);
    integrator->setKernels(kernels, solution, materials);

    PYLITH_METHOD_END;
} // _setKernelsJacobian


// ------------------------------------------------------------------------------------------------
// Find points with impulses.
void
pylith::faults::_FaultCohesiveImpulses::findImpulsePoints(int_array* impulsePoints,
                                                          const pylith::topology::Field& auxiliaryField,
                                                          const double threshold) {
    PYLITH_METHOD_BEGIN;
    // PYLITH_JOURNAL_DEBUG("_updateSlip(auxiliaryField="<<auxiliaryField<<", impulseStep="<<impulseStep<<")");

    PetscErrorCode err = 0;

    PetscSection auxiliaryFieldSection = auxiliaryField.getGlobalSection();assert(auxiliaryFieldSection);
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryFieldSection, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    pylith::topology::VecVisitorMesh auxiliaryVisitor(auxiliaryField, "slip");
    if (pStart == pEnd) {
        PYLITH_METHOD_END;
    } // if

    PylithScalar* auxiliaryArray = auxiliaryVisitor.localArray();assert(auxiliaryArray);

    // Use global section to get number of degrees of freedom so that we don't double count.
    PetscSection slipSectionGlobal = NULL;
    const int slipIndex = auxiliaryField.getSubfieldInfo("slip").index;
    err = PetscSectionGetField(auxiliaryFieldSection, slipIndex, &slipSectionGlobal);PYLITH_CHECK_ERROR(err);

    size_t count = 0;
    for (PetscInt point = pStart; point < pEnd; ++point) {
        PetscInt slipDof = 0;
        err = PetscSectionGetDof(slipSectionGlobal, point, &slipDof);PYLITH_CHECK_ERROR(err);

        const PetscInt slipOff = auxiliaryVisitor.sectionOffset(point);
        assert(auxiliaryVisitor.sectionDof(point) >= slipDof);
        for (PetscInt iDOF = 0; iDOF < slipDof; ++iDOF) {
            if (auxiliaryArray[slipOff+iDOF] >= threshold) {
                ++count;
                break;
            } // if
        } // for
    } // for

    impulsePoints->resize(count);
    size_t index = 0;
    for (PetscInt point = pStart; point < pEnd; ++point) {
        PetscInt slipDof = 0;
        err = PetscSectionGetDof(slipSectionGlobal, point, &slipDof);PYLITH_CHECK_ERROR(err);

        const PetscInt slipOff = auxiliaryVisitor.sectionOffset(point);
        for (PetscInt iDOF = 0; iDOF < slipDof; ++iDOF) {
            if (auxiliaryArray[slipOff+iDOF] >= threshold) {
                (*impulsePoints)[index++] = point;
                break;
            } // if
        } // for
    } // for

    PYLITH_METHOD_END;
} // findImpulsePoints


// End of file
