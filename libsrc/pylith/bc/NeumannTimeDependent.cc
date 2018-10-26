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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "NeumannTimeDependent.hh" // implementation of object methods

#include "TimeDependentAuxiliaryFactory.hh" // USES TimeDependentAuxiliaryFactory

#include "pylith/fekernels/NeumannTimeDependent.hh" // USES NeumannTimeDependent kernels

#include "pylith/feassemble/IntegratorBoundary.hh" // USES IntegratorBoundary
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorBoundary::ResidualKernels ResidualKernels;

// ---------------------------------------------------------------------------------------------------------------------
namespace pylith {
    namespace bc {
        class _NeumannTimeDependent {
            // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /** Set kernels for RHS residual.
             *
             * @param[out] integrator Integrator for boundary condition.
             * @param[in] bc Neumann time-dependent boundary condition.
             * @param[in] solution Solution field.
             */
            static
            void setKernelsRHSResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                       const pylith::bc::NeumannTimeDependent& bc,
                                       const pylith::topology::Field& solution);

        };

        // _NeumannTimeDependent

    } // bc
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::NeumannTimeDependent::NeumannTimeDependent(void) :
    _dbTimeHistory(NULL),
    _auxiliaryFactory(new pylith::bc::TimeDependentAuxiliaryFactory(pylith::bc::TimeDependentAuxiliaryFactory::TANGENTIAL_NORMAL)),
    _scaleName("pressure"),
    _useInitial(true),
    _useRate(false),
    _useTimeHistory(false) {
    PyreComponent::name("neumanntimedependent");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::NeumannTimeDependent::~NeumannTimeDependent(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::NeumannTimeDependent::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    BoundaryCondition::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
    _dbTimeHistory = NULL; // :KLUDGE: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set time history database.
void
pylith::bc::NeumannTimeDependent::setTimeHistoryDB(spatialdata::spatialdb::TimeHistory* th) {
    _dbTimeHistory = th;
} // setTimeHistoryDB


// ---------------------------------------------------------------------------------------------------------------------
// Get time history database.
const spatialdata::spatialdb::TimeHistory*
pylith::bc::NeumannTimeDependent::getTimeHistoryDB(void) {
    return _dbTimeHistory;
} // getTimeHistoryDB


// ---------------------------------------------------------------------------------------------------------------------
// Use initial value term in time history expression.
void
pylith::bc::NeumannTimeDependent::useInitial(const bool value) { // useInitial
    PYLITH_COMPONENT_DEBUG("useInitial(value="<<value<<")");

    _useInitial = value;
} // useInitial


// ---------------------------------------------------------------------------------------------------------------------
// Get flag associated with using initial value term in time history expression.
bool
pylith::bc::NeumannTimeDependent::useInitial(void) const { // useInitial
    return _useInitial;
} // useInitial


// ---------------------------------------------------------------------------------------------------------------------
// Use rate value term in time history expression.
void
pylith::bc::NeumannTimeDependent::useRate(const bool value) { // useRate
    PYLITH_COMPONENT_DEBUG("useRate(value="<<value<<")");

    _useRate = value;
} // useRate


// ---------------------------------------------------------------------------------------------------------------------
// Get flag associated with using rate value term in time history expression.
bool
pylith::bc::NeumannTimeDependent::useRate(void) const { // useRate
    return _useRate;
} // useRate


// ---------------------------------------------------------------------------------------------------------------------
// Use time history term in time history expression.
void
pylith::bc::NeumannTimeDependent::useTimeHistory(const bool value) { // useTimeHistory
    PYLITH_COMPONENT_DEBUG("useTimeHistory(value="<<value<<")");

    _useTimeHistory = value;
} // useTimeHistory


// ---------------------------------------------------------------------------------------------------------------------
// Get flag associated with using time history term in time history expression.
bool
pylith::bc::NeumannTimeDependent::useTimeHistory(void) const { // useTimeHistory
    return _useTimeHistory;
} // useTimeHistory


// ---------------------------------------------------------------------------------------------------------------------
// Name of scale associated with Neumann boundary condition (e.g., pressure for elasticity).
void
pylith::bc::NeumannTimeDependent::setScaleName(const char* value) {
    if (( value == std::string("length")) ||
        ( value == std::string("time")) ||
        ( value == std::string("pressure")) ||
        ( value == std::string("density")) ||
        ( value == std::string("pressure")) ) {
        _scaleName = value;
    } else {
        std::ostringstream msg;
        msg << "Unknown name of scale ("<<value<<") for Neumann boundary condition '" << _boundaryLabel << "'.";
        throw std::runtime_error(msg.str());
    } // if
} // setScaleName


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::bc::NeumannTimeDependent::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.label()<<")");

    pylith::feassemble::IntegratorBoundary* integrator = new pylith::feassemble::IntegratorBoundary(this);assert(integrator);
    integrator->setMarkerLabel(getMarkerLabel());

    _NeumannTimeDependent::setKernelsRHSResidual(integrator, *this, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
pylith::feassemble::Constraint*
pylith::bc::NeumannTimeDependent::createConstraint(const pylith::topology::Field& solution) {
    return NULL;
} // createConstraint


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::bc::NeumannTimeDependent::createAuxiliaryField(const pylith::topology::Field& solution,
                                                       const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.label()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->label("NeumannTimeDependent auxiliary field");

    assert(_auxiliaryFactory);
    assert(_normalizer);
    pylith::topology::Field::Description description = solution.subfieldInfo(_subfieldName.c_str()).description;
    if (_scaleName == std::string("pressure")) {
        description.scale = _normalizer->pressureScale();
    } else if (_scaleName == std::string("velocity")) {
        description.scale = _normalizer->lengthScale() / _normalizer->pressureScale();
    } else if (_scaleName == std::string("length")) {
        description.scale = _normalizer->pressureScale();
    } else if (_scaleName == std::string("time")) {
        description.scale = _normalizer->pressureScale();
    } else if (_scaleName == std::string("debsuty")) {
        description.scale = _normalizer->pressureScale();
    } else {
        std::ostringstream msg;
        msg << "Unknown name of scale ("<<_scaleName<<") for Neumann boundary condition for '" << _boundaryLabel << "'.";
        PYLITH_COMPONENT_ERROR(msg.str());
        throw std::logic_error(msg.str());
    } // if/else
    _auxiliaryFactory->initialize(auxiliaryField, *_normalizer, solution.spaceDim(), &description);

    // :ATTENTION: The order of the factory methods must match the order of the auxiliary subfields in the FE kernels.

    if (_useInitial) {
        _auxiliaryFactory->addInitialAmplitude();
    } // if
    if (_useRate) {
        _auxiliaryFactory->addRateAmplitude();
        _auxiliaryFactory->addRateStartTime();
    } // _useRate
    if (_useTimeHistory) {
        _auxiliaryFactory->addTimeHistoryAmplitude();
        _auxiliaryFactory->addTimeHistoryStartTime();
        _auxiliaryFactory->addTimeHistoryValue();
    } // _useTimeHistory

    auxiliaryField->subfieldsSetup();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->zeroLocal();

    assert(_auxiliaryFactory);
    _auxiliaryFactory->setValuesFromDB();

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::bc::NeumannTimeDependent::createDerivedField(const pylith::topology::Field& solution,
                                                     const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.label()<<", domainMesh=)"<<typeid(domainMesh).name()<<") empty method");

    PYLITH_METHOD_RETURN(NULL);
} // createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary subfields at beginning of time step.
void
pylith::bc::NeumannTimeDependent::prestep(pylith::topology::Field* auxiliaryField,
                                          const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("prestep(auxiliaryField="<<auxiliaryField<<", t="<<t<<")");

    if (_useTimeHistory) {
        assert(_normalizer);
        const PylithScalar timeScale = _normalizer->timeScale();
        TimeDependentAuxiliaryFactory::updateAuxiliaryField(auxiliaryField, t, timeScale, _dbTimeHistory);
    } // if

    PYLITH_METHOD_END;
} // prestep


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::bc::NeumannTimeDependent::_getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // _getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::bc::NeumannTimeDependent::_updateKernelConstants(const PylithReal dt) {
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
// Set kernels for RHS residual G(t,s).
void
pylith::bc::_NeumannTimeDependent::setKernelsRHSResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                                         const pylith::bc::NeumannTimeDependent& bc,
                                                         const topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    //PYLITH_COMPONENT_DEBUG("setKernelsRHSResidual(integrator="<<integrator<<", solution="<<solution.label()<<")");

    PetscBdPointFunc g0 = NULL;
    PetscBdPointFunc g1 = NULL;

    const pylith::topology::Field::VectorFieldEnum fieldType = solution.subfieldInfo(bc.getSubfieldName()).description.vectorFieldType;
    const bool isScalarField = fieldType == pylith::topology::Field::SCALAR;

    const int bitInitial = bc.useInitial() ? 0x1 : 0x0;
    const int bitRate = bc.useRate() ? 0x2 : 0x0;
    const int bitTimeHistory = bc.useTimeHistory() ? 0x4 : 0x0;
    const int bitUse = bitInitial | bitRate | bitTimeHistory;
    switch (bitUse) {
    case 0x1:
        g0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::g0_initial_scalar :
             pylith::fekernels::NeumannTimeDependent::g0_initial_vector;
        break;
    case 0x2:
        g0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::g0_rate_scalar :
             pylith::fekernels::NeumannTimeDependent::g0_rate_vector;
        break;
    case 0x4:
        g0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::g0_timeHistory_scalar :
             pylith::fekernels::NeumannTimeDependent::g0_timeHistory_vector;
        break;
    case 0x3:
        g0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::g0_initialRate_scalar :
             pylith::fekernels::NeumannTimeDependent::g0_initialRate_scalar;
        break;
    case 0x5:
        g0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::g0_initialTimeHistory_scalar :
             pylith::fekernels::NeumannTimeDependent::g0_initialTimeHistory_scalar;
        break;
    case 0x6:
        g0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::g0_rateTimeHistory_scalar :
             pylith::fekernels::NeumannTimeDependent::g0_rateTimeHistory_vector;
        break;
    case 0x7:
        g0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::g0_initialRateTimeHistory_scalar :
             pylith::fekernels::NeumannTimeDependent::g0_initialRateTimeHistory_vector;
        break;
    case 0x0:
        //PYLITH_COMPONENT_WARNING("Neumann time-dependent BC provides no values.");
        break;
    default:
        throw std::logic_error("Unknown combination of flags for Neumann time-dependent BC terms.");
    } // switch

    std::vector<ResidualKernels> kernels(1);
    kernels[0] = ResidualKernels(bc.getSubfieldName(), g0, g1);

    assert(integrator);
    integrator->setKernelsRHSResidual(kernels);

    PYLITH_METHOD_END;
} // setKernelsRHSResidual


// End of file
