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

#include "FaultCohesiveKinPoro.hh" // implementation of object methods

#include "pylith/faults/KinSrc.hh"                    // USES KinSrc
#include "pylith/faults/AuxiliaryFactoryKinematic.hh" // USES AuxiliaryFactoryKinematic
#include "pylith/feassemble/IntegratorInterface.hh"   // USES IntegratorInterface
#include "pylith/feassemble/InterfacePatches.hh"      // USES InterfacePatches
#include "pylith/feassemble/ConstraintSimple.hh"      // USES ConstraintSimple

#include "pylith/topology/Mesh.hh"        // USES Mesh
#include "pylith/topology/Field.hh"       // USES Field
#include "pylith/topology/FieldOps.hh"    // USES FieldOps::checkDiscretization()
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/fekernels/FaultCohesiveKinPoro.hh" // USES FaultCohesiveKinPoro

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh"    // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh"   // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensionalizer
#include "spatialdata/spatialdb/SpatialDB.hh"  // USES SpatialDB

#include <cmath>     // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring>   // USES strlen()
#include <cstdlib>   // USES atoi()
#include <cassert>   // USES assert()
#include <sstream>   // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error
#include <typeinfo>  // USES typeid()

// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorInterface::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorInterface::JacobianKernels JacobianKernels;

// ------------------------------------------------------------------------------------------------
namespace pylith
{
    namespace faults
    {
        class _FaultCohesiveKinPoro {
            // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////
        public:
            static const char *pyreComponent;
        };
        const char* _FaultCohesiveKinPoro::pyreComponent = "faultcohesivekinporo";

        // _FaultCohesiveKinPoro

    } // faults
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveKinPoro::FaultCohesiveKinPoro(void) : _auxiliaryFactory(new pylith::faults::AuxiliaryFactoryKinematic),
                                                                   _slipVecRupture(NULL),
                                                                   _slipVecTotal(NULL)
{
    pylith::utils::PyreComponent::setName(_FaultCohesiveKinPoro::pyreComponent);
} // constructor

// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveKinPoro::~FaultCohesiveKinPoro(void)
{
    deallocate();
} // destructor

// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void pylith::faults::FaultCohesiveKinPoro::deallocate(void)
{
    FaultCohesive::deallocate();

    PetscErrorCode err = VecDestroy(&_slipVecRupture);
    PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_slipVecTotal);
    PYLITH_CHECK_ERROR(err);
    delete _auxiliaryFactory;
    _auxiliaryFactory = NULL;
    _ruptures.clear(); // :TODO: Use shared pointers for earthquake ruptures
} // deallocate

// ------------------------------------------------------------------------------------------------
// Set kinematic earthquake ruptures.
void pylith::faults::FaultCohesiveKinPoro::setEqRuptures(const char *const *names,
                                                         const int numNames,
                                                         KinSrc **ruptures,
                                                         const int numRuptures)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setEqRuptures(names=" << names << ", numNames=" << numNames << ", ruptures=" << ruptures << ", numRuptures=" << numRuptures << ")");

    assert(numNames == numRuptures);

    // :TODO: Use shared pointers for earthquake ruptures
    _ruptures.clear();
    for (int i = 0; i < numRuptures; ++i)
    {
        if (!ruptures[i])
        {
            std::ostringstream msg;
            msg << "Null earthquake rupture object for earthquake rupture '" << names[i] << "'.";
            throw std::runtime_error(msg.str());
        } // if
        _ruptures[std::string(names[i])] = ruptures[i];
    } // for

    PYLITH_METHOD_END;
} // setEqRuptures

// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void pylith::faults::FaultCohesiveKinPoro::verifyConfiguration(const pylith::topology::Field &solution) const
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution=" << solution.getLabel() << ")");

    if (!solution.hasSubfield("lagrange_multiplier_fault"))
    {
        std::ostringstream msg;
        msg << "Cannot find 'lagrange_multiplier_fault' subfield in solution field for fault implementation in component '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    if (!solution.hasSubfield("mu_multiplier_fault"))
    {
        std::ostringstream msg;
        msg << "Cannot find 'mu_multiplier_fault' subfield in solution field for fault implementation in component '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if    

    switch (_formulation)
    {
    case QUASISTATIC:
        if (!solution.hasSubfield("displacement"))
        {
            std::ostringstream msg;
            msg << "Cannot find 'displacement' subfield in solution field for fault implementation in component '"
                << PyreComponent::getIdentifier() << "'.";
            throw std::runtime_error(msg.str());
        } // if
        break;
    case DYNAMIC_IMEX:
        if (!solution.hasSubfield("velocity"))
        {
            std::ostringstream msg;
            msg << "Cannot find 'velocity' subfield in solution field for fault implementation in component '"
                << PyreComponent::getIdentifier() << "'.";
            throw std::runtime_error(msg.str());
        } // if
        break;
    case DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' formulation. Use 'dynamic_imex'.");
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    PYLITH_METHOD_END;
} // verifyConfiguration

// ------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator *
pylith::faults::FaultCohesiveKinPoro::createIntegrator(const pylith::topology::Field &solution)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution=" << solution.getLabel() << ")");

    pylith::feassemble::IntegratorInterface *integrator = new pylith::feassemble::IntegratorInterface(this);
    assert(integrator);
    integrator->setLabelValue(getInterfaceId());
    integrator->setSurfaceMarkerLabel(getSurfaceMarkerLabel());

    pylith::feassemble::InterfacePatches *patches =
        pylith::feassemble::InterfacePatches::createMaterialPairs(this, solution.getDM());
    integrator->setIntegrationPatches(patches);

    _setKernelsResidual(integrator, solution);
    _setKernelsJacobian(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator

// ------------------------------------------------------------------------------------------------
// Create constraint for buried fault edges and faces.
std::vector<pylith::feassemble::Constraint*>
pylith::faults::FaultCohesiveKinPoro::createConstraints(const pylith::topology::Field &solution)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createConstraints(solution=" << solution.getLabel() << ")");

    std::vector<pylith::feassemble::Constraint*> constraintArray;

    if (0 == strlen(getBuriedEdgesMarkerLabel())) {
        PYLITH_METHOD_RETURN(constraintArray);
    } // if

    // Lagrange Multipliers
    const char *lagrangeName = "lagrange_multiplier_fault";
    const PylithInt numComponents = solution.getSpaceDim();

    pylith::int_array constrainedDOFLagrange;
    constrainedDOFLagrange.resize(numComponents);
    for (int c = 0; c < numComponents; ++c)
    {
        constrainedDOFLagrange[c] = c;
    }
    // Make new label for cohesive edges and faces
    PetscDM dm = solution.getDM();
    PetscDMLabel buriedLabel = NULL;
    PetscDMLabel buriedCohesiveLabel = NULL;
    PetscIS pointIS = NULL;
    const PetscInt *points = NULL;
    PetscInt n;
    std::ostringstream labelstream;
    labelstream << getBuriedEdgesMarkerLabel() << "_cohesive";
    std::string labelname = labelstream.str();
    PetscErrorCode err;

    err = DMCreateLabel(dm, labelname.c_str());PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dm, getBuriedEdgesMarkerLabel(), &buriedLabel);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dm, labelname.c_str(), &buriedCohesiveLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(buriedLabel, 1, &pointIS);PYLITH_CHECK_ERROR(err);
    err = ISGetLocalSize(pointIS, &n);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    for (int p = 0; p < n; ++p) {
        const PetscInt *support = NULL;
        PetscInt supportSize;

        err = DMPlexGetSupportSize(dm, points[p], &supportSize);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetSupport(dm, points[p], &support);PYLITH_CHECK_ERROR(err);
        for (int s = 0; s < supportSize; ++s) {
            DMPolytopeType ct;
            const PetscInt spoint = support[s];

            err = DMPlexGetCellType(dm, spoint, &ct);
            PYLITH_CHECK_ERROR(err);
            if ((ct == DM_POLYTOPE_SEG_PRISM_TENSOR) || (ct == DM_POLYTOPE_POINT_PRISM_TENSOR)) {
                const PetscInt *cone = NULL;
                PetscInt coneSize;

                err = DMPlexGetConeSize(dm, spoint, &coneSize);PYLITH_CHECK_ERROR(err);
                err = DMPlexGetCone(dm, spoint, &cone);PYLITH_CHECK_ERROR(err);
                for (int c = 0; c < coneSize; ++c) {
                    PetscInt val;
                    err = DMLabelGetValue(buriedLabel, cone[c], &val);PYLITH_CHECK_ERROR(err);
                    if (val >= 0) {
                        err = DMLabelSetValue(buriedCohesiveLabel, spoint, 1);PYLITH_CHECK_ERROR(err);
                        break;
                    }
                } // for
            } // if
        } // for
    } // for

    pylith::feassemble::ConstraintSimple *constraintLagrange = new pylith::feassemble::ConstraintSimple(this);assert(constraintLagrange);
    constraintLagrange->setMarkerLabel(labelname.c_str());
    err = PetscObjectViewFromOptions((PetscObject)buriedLabel, NULL, "-buried_edge_label_view");
    err = PetscObjectViewFromOptions((PetscObject)buriedCohesiveLabel, NULL, "-buried_cohesive_edge_label_view");
    constraintLagrange->setConstrainedDOF(&constrainedDOFLagrange[0], constrainedDOFLagrange.size());
    constraintLagrange->setSubfieldName(lagrangeName);
    constraintLagrange->setUserFn(_zero);

    // "Mu" multipliers (test case same as lagrange)
    const char *muName = "mu_multiplier_fault";

    pylith::int_array constrainedDOFMu;
    constrainedDOFMu.resize(numComponents);
    for (int c = numComponents; c < numComponents*2; ++c)
    {
        constrainedDOFMu[c] = c;
    }
     // Make new label for cohesive edges and faces
    //PetscDM dm = solution.getDM();
    PetscDMLabel buriedLabelMu = NULL;
    PetscDMLabel buriedCohesiveLabelMu = NULL;
    PetscIS pointISMu = NULL;
    const PetscInt *pointsMu = NULL;
    PetscInt nMu;
    std::ostringstream labelstreamMu;
    labelstreamMu << getBuriedEdgesMarkerLabel() << "_cohesive";
    std::string labelnameMu = labelstreamMu.str();
    //PetscErrorCode err;

    err = DMCreateLabel(dm, labelnameMu.c_str());PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dm, getBuriedEdgesMarkerLabel(), &buriedLabelMu);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dm, labelnameMu.c_str(), &buriedCohesiveLabelMu);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(buriedLabelMu, 1, &pointISMu);PYLITH_CHECK_ERROR(err);
    err = ISGetLocalSize(pointISMu, &nMu);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(pointISMu, &pointsMu);PYLITH_CHECK_ERROR(err);
    for (int p = 0; p < nMu; ++p) {
        const PetscInt *supportMu = NULL;
        PetscInt supportSizeMu;

        err = DMPlexGetSupportSize(dm, pointsMu[p], &supportSizeMu);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetSupport(dm, pointsMu[p], &supportMu);PYLITH_CHECK_ERROR(err);
        for (int s = 0; s < supportSizeMu; ++s) {
            DMPolytopeType ctMu;
            const PetscInt spointMu = supportMu[s];

            err = DMPlexGetCellType(dm, spointMu, &ctMu);
            PYLITH_CHECK_ERROR(err);
            if ((ctMu == DM_POLYTOPE_SEG_PRISM_TENSOR) || (ctMu == DM_POLYTOPE_POINT_PRISM_TENSOR)) {
                const PetscInt *coneMu = NULL;
                PetscInt coneSizeMu;

                err = DMPlexGetConeSize(dm, spointMu, &coneSizeMu);PYLITH_CHECK_ERROR(err);
                err = DMPlexGetCone(dm, spointMu, &coneMu);PYLITH_CHECK_ERROR(err);
                for (int c = 0; c < coneSizeMu; ++c) {
                    PetscInt valMu;
                    err = DMLabelGetValue(buriedLabelMu, coneMu[c], &valMu);PYLITH_CHECK_ERROR(err);
                    if (valMu >= 0) {
                        err = DMLabelSetValue(buriedCohesiveLabelMu, spointMu, 1);PYLITH_CHECK_ERROR(err);
                        break;
                    }
                } // for
            } // if
        } // for
    } // for

    pylith::feassemble::ConstraintSimple *constraintMu = new pylith::feassemble::ConstraintSimple(this);assert(constraintMu);
    constraintMu->setMarkerLabel(labelnameMu.c_str());
    err = PetscObjectViewFromOptions((PetscObject)buriedLabelMu, NULL, "-buried_edge_label_view");
    err = PetscObjectViewFromOptions((PetscObject)buriedCohesiveLabelMu, NULL, "-buried_cohesive_edge_label_view");
    constraintMu->setConstrainedDOF(&constrainedDOFMu[0], constrainedDOFMu.size());
    constraintMu->setSubfieldName(muName);
    constraintMu->setUserFn(_zero);

    // Package constraints and exit
    constraintArray.resize(2);
    constraintArray[0] = constraintLagrange;
    constraintArray[1] = constraintMu;    
    PYLITH_METHOD_RETURN(constraintArray);
} // createConstraints

// ------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field *
pylith::faults::FaultCohesiveKinPoro::createAuxiliaryField(const pylith::topology::Field &solution,
                                                           const pylith::topology::Mesh &domainMesh)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution=" << solution.getLabel() << ", domainMesh=)" << typeid(domainMesh).name() << ")");

    assert(_normalizer);

    pylith::topology::Field *auxiliaryField = new pylith::topology::Field(domainMesh);
    assert(auxiliaryField);
    auxiliaryField->setLabel("FaultCohesiveKinPoro auxiliary field");

    // Set default discretization of auxiliary subfields to match lagrange_multiplier_fault subfield in solution.
    assert(_auxiliaryFactory);
    const pylith::topology::FieldBase::Discretization &discretization = solution.getSubfieldInfo("lagrange_multiplier_fault").fe;
    const PylithInt dimension = -1;
    const bool isFaultOnly = false;
    _auxiliaryFactory->setSubfieldDiscretization("default", discretization.basisOrder, discretization.quadOrder, dimension,
                                                 isFaultOnly, discretization.cellBasis, discretization.feSpace,
                                                 discretization.isBasisContinuous);

    assert(_auxiliaryFactory);
    assert(_normalizer);
    _auxiliaryFactory->initialize(auxiliaryField, *_normalizer, solution.getSpaceDim());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    switch (_formulation)
    {
    case QUASISTATIC:
        _auxiliaryFactory->addSlip(); // 0`
        break;
    case DYNAMIC_IMEX:
        _auxiliaryFactory->addSlipAcceleration(); // 0
        break;
    case DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' time-stepping formulation. Use 'dynamic_imex'.");
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    // We don't populate the auxiliary field via a spatial database, because they will be set from the earthquake
    // rupture.

    // Initialize auxiliary fields for kinematic ruptures.
    assert(auxiliaryField);
    const srcs_type::const_iterator rupturesEnd = _ruptures.end();
    for (srcs_type::iterator r_iter = _ruptures.begin(); r_iter != rupturesEnd; ++r_iter)
    {
        KinSrc *src = r_iter->second;
        assert(src);
        src->initialize(*auxiliaryField, *_normalizer, solution.getMesh().getCoordSys());
    } // for

    // Create local PETSc vector to hold current slip.
    PetscErrorCode err = 0;
    err = DMCreateLocalVector(auxiliaryField->getDM(), &_slipVecRupture);
    PYLITH_CHECK_ERROR(err);
    err = DMCreateLocalVector(auxiliaryField->getDM(), &_slipVecTotal);
    PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField

// ------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field *
pylith::faults::FaultCohesiveKinPoro::createDerivedField(const pylith::topology::Field &solution,
                                                         const pylith::topology::Mesh &domainMesh)
{
    return NULL;
} // createDerivedField

// ------------------------------------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void pylith::faults::FaultCohesiveKinPoro::updateAuxiliaryField(pylith::topology::Field *auxiliaryField,
                                                                const double t)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateAuxiliaryField(auxiliaryField=" << auxiliaryField << ", t=" << t << ")");

    assert(auxiliaryField);
    assert(_normalizer);

    switch (_formulation)
    {
    case QUASISTATIC:
        this->_updateSlip(auxiliaryField, t);
        break;
    case DYNAMIC_IMEX:
        this->_updateSlipAcceleration(auxiliaryField, t);
        break;
    case DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' formulation. Use 'dynamic_imex'.");
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    PYLITH_METHOD_END;
} // updateAuxiliaryField

// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory *
pylith::faults::FaultCohesiveKinPoro::_getAuxiliaryFactory(void)
{
    return _auxiliaryFactory;
} // _getAuxiliaryFactory

// ------------------------------------------------------------------------------------------------
// Update kernel constants.
void pylith::faults::FaultCohesiveKinPoro::_updateKernelConstants(const PylithReal dt)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelConstants(dt=" << dt << ")");

    if (6 != _kernelConstants.size())
    {
        _kernelConstants.resize(6);
    }
    _kernelConstants[0] = _refDir1[0];
    _kernelConstants[1] = _refDir1[1];
    _kernelConstants[2] = _refDir1[2];
    _kernelConstants[3] = _refDir2[0];
    _kernelConstants[4] = _refDir2[1];
    _kernelConstants[5] = _refDir2[2];

    PYLITH_METHOD_END;
} // _updateKernelConstants

// ------------------------------------------------------------------------------------------------
// Update slip subfield in auxiliary field at beginning of time step.
void pylith::faults::FaultCohesiveKinPoro::_updateSlip(pylith::topology::Field *auxiliaryField,
                                                       const double t)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateSlip(auxiliaryField=" << auxiliaryField << ", t=" << t << ")");

    assert(auxiliaryField);
    assert(_normalizer);

    // Update slip subfield at current time step
    PetscErrorCode err = VecSet(_slipVecTotal, 0.0);
    PYLITH_CHECK_ERROR(err);
    const srcs_type::const_iterator rupturesEnd = _ruptures.end();
    for (srcs_type::iterator r_iter = _ruptures.begin(); r_iter != rupturesEnd; ++r_iter)
    {
        err = VecSet(_slipVecRupture, 0.0);
        PYLITH_CHECK_ERROR(err);

        KinSrc *src = r_iter->second;
        assert(src);
        src->updateSlip(_slipVecRupture, auxiliaryField, t, _normalizer->getTimeScale());
        err = VecAYPX(_slipVecTotal, 1.0, _slipVecRupture);
    } // for

    // Transfer slip values from local PETSc slip vector to fault auxiliary field.
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryField->getLocalSection(), &pStart, &pEnd);
    PYLITH_CHECK_ERROR(err);

    pylith::topology::VecVisitorMesh auxiliaryVisitor(*auxiliaryField, "slip");
    PylithScalar *auxiliaryArray = auxiliaryVisitor.localArray();

    const PylithScalar *slipArray = NULL;
    err = VecGetArrayRead(_slipVecTotal, &slipArray);
    PYLITH_CHECK_ERROR(err);

    for (PetscInt p = pStart, iSlip = 0; p < pEnd; ++p)
    {
        const PetscInt slipDof = auxiliaryVisitor.sectionDof(p);
        const PetscInt slipOff = auxiliaryVisitor.sectionOffset(p);
        for (PetscInt iDof = 0; iDof < slipDof; ++iDof, ++iSlip)
        {
            auxiliaryArray[slipOff + iDof] = slipArray[iSlip];
        } // for
    }     // for
    err = VecRestoreArrayRead(_slipVecTotal, &slipArray);
    PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state())
    {
        auxiliaryField->view("Fault auxiliary field after setting slip.");
    } // if

    PYLITH_METHOD_END;
} // _updateSlip

// ------------------------------------------------------------------------------------------------
// Update slip rate subfield in auxiliary field at beginning of time step.
void pylith::faults::FaultCohesiveKinPoro::_updateSlipRate(pylith::topology::Field *auxiliaryField,
                                                           const double t)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_updateSlipRate(auxiliaryField=" << auxiliaryField << ", t=" << t << ")");

    assert(auxiliaryField);
    assert(_normalizer);

    // Update slip rate subfield at current time step
    PetscErrorCode err = VecSet(_slipVecTotal, 0.0);
    PYLITH_CHECK_ERROR(err);
    const srcs_type::const_iterator rupturesEnd = _ruptures.end();
    for (srcs_type::iterator r_iter = _ruptures.begin(); r_iter != rupturesEnd; ++r_iter)
    {
        err = VecSet(_slipVecRupture, 0.0);
        PYLITH_CHECK_ERROR(err);

        KinSrc *src = r_iter->second;
        assert(src);
        src->updateSlipRate(_slipVecRupture, auxiliaryField, t, _normalizer->getTimeScale());
        err = VecAYPX(_slipVecTotal, 1.0, _slipVecRupture);
    } // for

    // Transfer slip values from local PETSc slip vector to fault auxiliary field.
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryField->getLocalSection(), &pStart, &pEnd);
    PYLITH_CHECK_ERROR(err);

    pylith::topology::VecVisitorMesh auxiliaryVisitor(*auxiliaryField, "slip_rate");
    PylithScalar *auxiliaryArray = auxiliaryVisitor.localArray();

    const PylithScalar *slipRateArray = NULL;
    err = VecGetArrayRead(_slipVecTotal, &slipRateArray);
    PYLITH_CHECK_ERROR(err);

    for (PetscInt p = pStart, iSlipRate = 0; p < pEnd; ++p)
    {
        const PetscInt slipRateDof = auxiliaryVisitor.sectionDof(p);
        const PetscInt slipRateOff = auxiliaryVisitor.sectionOffset(p);
        for (PetscInt iDof = 0; iDof < slipRateDof; ++iDof, ++iSlipRate)
        {
            auxiliaryArray[slipRateOff + iDof] = slipRateArray[iSlipRate];
        } // for
    }     // for
    err = VecRestoreArrayRead(_slipVecTotal, &slipRateArray);
    PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state())
    {
        auxiliaryField->view("Fault auxiliary field after setting slip rate.");
    } // if

    PYLITH_METHOD_END;
} // _updateSlipRate

// ------------------------------------------------------------------------------------------------
// Update slip acceleration subfield in auxiliary field at beginning of time step.
void pylith::faults::FaultCohesiveKinPoro::_updateSlipAcceleration(pylith::topology::Field *auxiliaryField,
                                                                   const double t)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_updateSlipAcceleration(auxiliaryField=" << auxiliaryField << ", t=" << t << ")");

    assert(auxiliaryField);
    assert(_normalizer);

    // Update slip acceleration subfield at current time step
    PetscErrorCode err = VecSet(_slipVecTotal, 0.0);
    PYLITH_CHECK_ERROR(err);
    const srcs_type::const_iterator rupturesEnd = _ruptures.end();
    for (srcs_type::iterator r_iter = _ruptures.begin(); r_iter != rupturesEnd; ++r_iter)
    {
        err = VecSet(_slipVecRupture, 0.0);
        PYLITH_CHECK_ERROR(err);

        KinSrc *src = r_iter->second;
        assert(src);
        src->updateSlipAcc(_slipVecRupture, auxiliaryField, t, _normalizer->getTimeScale());
        err = VecAYPX(_slipVecTotal, 1.0, _slipVecRupture);
    } // for

    // Transfer slip values from local PETSc slip vector to fault auxiliary field.
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryField->getLocalSection(), &pStart, &pEnd);
    PYLITH_CHECK_ERROR(err);

    pylith::topology::VecVisitorMesh auxiliaryVisitor(*auxiliaryField, "slip_acceleration");
    PylithScalar *auxiliaryArray = auxiliaryVisitor.localArray();

    const PylithScalar *slipAccArray = NULL;
    err = VecGetArrayRead(_slipVecTotal, &slipAccArray);
    PYLITH_CHECK_ERROR(err);

    for (PetscInt p = pStart, iSlipAcc = 0; p < pEnd; ++p)
    {
        const PetscInt slipAccDof = auxiliaryVisitor.sectionDof(p);
        const PetscInt slipAccOff = auxiliaryVisitor.sectionOffset(p);
        for (PetscInt iDof = 0; iDof < slipAccDof; ++iDof, ++iSlipAcc)
        {
            auxiliaryArray[slipAccOff + iDof] = slipAccArray[iSlipAcc];
        } // for
    } // for
    err = VecRestoreArrayRead(_slipVecTotal, &slipAccArray);
    PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state())
    {
        auxiliaryField->view("Fault auxiliary field after setting slip acceleration.");
    } // if

    PYLITH_METHOD_END;
} // _updateSlipAcceleration

// ------------------------------------------------------------------------------------------------
// Set kernels for residual.
void pylith::faults::FaultCohesiveKinPoro::_setKernelsResidual(pylith::feassemble::IntegratorInterface *integrator,
                                                               const pylith::topology::Field &solution) const
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsResidual(integrator=" << integrator << ", solution=" << solution.getLabel() << ")");
    typedef pylith::feassemble::IntegratorInterface integrator_t;

    std::vector<ResidualKernels> kernels;
    switch (_formulation)
    {
    case pylith::problems::Physics::QUASISTATIC:
    {
        // Elasticity equation (displacement) for negative side of the fault.
        const PetscBdPointFunc f0u_neg = pylith::fekernels::FaultCohesiveKinPoro::f0u_neg;
        const PetscBdPointFunc f1u_neg = NULL;

        // Elasticity equation (displacement) for positive side of the fault.
        const PetscBdPointFunc f0u_pos = pylith::fekernels::FaultCohesiveKinPoro::f0u_pos;
        const PetscBdPointFunc f1u_pos = NULL;

        // Fault slip constraint equation, first.
        const PetscBdPointFunc f0l = pylith::fekernels::FaultCohesiveKinPoro::f0l_u;
        const PetscBdPointFunc f1l = NULL;

        // Fault slip constraint equation, second.
        const PetscBdPointFunc f0mu = pylith::fekernels::FaultCohesiveKinPoro::f0mu_u;
        const PetscBdPointFunc f1mu = NULL;

        kernels.resize(4);
        kernels[0] = ResidualKernels("displacement", integrator_t::RESIDUAL_LHS, integrator_t::NEGATIVE_FACE,
                                     f0u_neg, f1u_neg);
        kernels[1] = ResidualKernels("displacement", integrator_t::RESIDUAL_LHS, integrator_t::POSITIVE_FACE,
                                     f0u_pos, f1u_pos);
        kernels[2] = ResidualKernels("lagrange_multiplier_fault", integrator_t::RESIDUAL_LHS, integrator_t::FAULT_FACE,
                                     f0l, f1l);
        kernels[3] = ResidualKernels("mu_multiplier_fault", integrator_t::RESIDUAL_LHS, integrator_t::FAULT_FACE,
                                     f0mu, f1mu);

        break;
    } // QUASISTATIC
    case pylith::problems::Physics::DYNAMIC_IMEX:
    {
        break;
    } // DYNAMIC_IMEX
    case pylith::problems::Physics::DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' formulation. Use 'dynamic_imex'.");
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation '" << _formulation << "'.");
    } // switch

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsResidual

// ------------------------------------------------------------------------------------------------
// Set kernels for Jacobian.
void pylith::faults::FaultCohesiveKinPoro::_setKernelsJacobian(pylith::feassemble::IntegratorInterface *integrator,
                                                               const pylith::topology::Field &solution) const
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsJacobian(integrator=" << integrator << ", solution=" << solution.getLabel() << ")");
    typedef pylith::feassemble::IntegratorInterface integrator_t;

    std::vector<JacobianKernels> kernels;
    switch (_formulation)
    {
    case QUASISTATIC:
    {
        const PetscBdPointJac Jf0ul_neg = pylith::fekernels::FaultCohesiveKinPoro::Jf0ul_neg;
        const PetscBdPointJac Jf1ul_neg = NULL;
        const PetscBdPointJac Jf2ul_neg = NULL;
        const PetscBdPointJac Jf3ul_neg = NULL;

        const PetscBdPointJac Jf0ul_pos = pylith::fekernels::FaultCohesiveKinPoro::Jf0ul_pos;
        const PetscBdPointJac Jf1ul_pos = NULL;
        const PetscBdPointJac Jf2ul_pos = NULL;
        const PetscBdPointJac Jf3ul_pos = NULL;

        const PetscBdPointJac Jf0umu_neg = pylith::fekernels::FaultCohesiveKinPoro::Jf0umu_neg;
        const PetscBdPointJac Jf1umu_neg = NULL;
        const PetscBdPointJac Jf2umu_neg = NULL;
        const PetscBdPointJac Jf3umu_neg = NULL;

        const PetscBdPointJac Jf0umu_pos = pylith::fekernels::FaultCohesiveKinPoro::Jf0umu_pos;
        const PetscBdPointJac Jf1umu_pos = NULL;
        const PetscBdPointJac Jf2umu_pos = NULL;
        const PetscBdPointJac Jf3umu_pos = NULL;

        const PetscBdPointJac Jf0lu = pylith::fekernels::FaultCohesiveKinPoro::Jf0lu;
        const PetscBdPointJac Jf1lu = NULL;
        const PetscBdPointJac Jf2lu = NULL;
        const PetscBdPointJac Jf3lu = NULL;

        const PetscBdPointJac Jf0muu = pylith::fekernels::FaultCohesiveKinPoro::Jf0muu;
        const PetscBdPointJac Jf1muu = NULL;
        const PetscBdPointJac Jf2muu = NULL;
        const PetscBdPointJac Jf3muu = NULL;

        kernels.resize(6);
        const char *nameDisplacement = "displacement";
        const char *nameLagrangeMultiplier = "lagrange_multiplier_fault";
        const char *nameMu = "mu_multiplier_fault";
        kernels[0] = JacobianKernels(nameDisplacement, nameLagrangeMultiplier, integrator_t::JACOBIAN_LHS,
                                     integrator_t::NEGATIVE_FACE, Jf0ul_neg, Jf1ul_neg, Jf2ul_neg, Jf3ul_neg);
        kernels[1] = JacobianKernels(nameDisplacement, nameLagrangeMultiplier, integrator_t::JACOBIAN_LHS,
                                     integrator_t::POSITIVE_FACE, Jf0ul_pos, Jf1ul_pos, Jf2ul_pos, Jf3ul_pos);
        kernels[2] = JacobianKernels(nameDisplacement, nameMu, integrator_t::JACOBIAN_LHS,
                                     integrator_t::NEGATIVE_FACE, Jf0umu_neg, Jf1umu_neg, Jf2umu_neg, Jf3umu_neg);
        kernels[3] = JacobianKernels(nameDisplacement, nameMu, integrator_t::JACOBIAN_LHS,
                                     integrator_t::POSITIVE_FACE, Jf0umu_pos, Jf1umu_pos, Jf2umu_pos, Jf3umu_pos);
        kernels[4] = JacobianKernels(nameLagrangeMultiplier, nameDisplacement, integrator_t::JACOBIAN_LHS,
                                     integrator_t::FAULT_FACE, Jf0lu, Jf1lu, Jf2lu, Jf3lu);
        kernels[5] = JacobianKernels(nameMu, nameDisplacement, integrator_t::JACOBIAN_LHS,
                                     integrator_t::FAULT_FACE, Jf0muu, Jf1muu, Jf2muu, Jf3muu);
        break;
    } // QUASISTATIC
    case pylith::problems::Physics::DYNAMIC_IMEX:
    {
        break;
    } // DYNAMIC_IMEX
    case pylith::problems::Physics::DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' formulation. Use 'dynamic_imex'.");
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation '" << _formulation << "'.");
    } // switch

    assert(integrator);
    integrator->setKernelsJacobian(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsJacobian

// End of file