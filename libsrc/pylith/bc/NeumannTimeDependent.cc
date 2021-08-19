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

#include "NeumannTimeDependent.hh" // implementation of object methods

#include "TimeDependentAuxiliaryFactory.hh" // USES TimeDependentAuxiliaryFactory

#include "pylith/fekernels/NeumannTimeDependent.hh" // USES NeumannTimeDependent kernels

#include "pylith/feassemble/IntegratorBoundary.hh" // USES IntegratorBoundary
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
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
             * @param[in] formulation Formulation for equations.
             */
            static
            void setKernelsResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                    const pylith::bc::NeumannTimeDependent& bc,
                                    const pylith::topology::Field& solution,
                                    const pylith::problems::Physics::FormulationEnum formulation);

            static const char* pyreComponent;

        }; // _NeumannTimeDependent
        const char* _NeumannTimeDependent::pyreComponent = "neumanntimedependent";

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
    PyreComponent::setName(_NeumannTimeDependent::pyreComponent);
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
    PYLITH_COMPONENT_DEBUG("setTimeHistoryDB(th"<<th<<")");

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
pylith::bc::NeumannTimeDependent::useInitial(const bool value) {
    PYLITH_COMPONENT_DEBUG("useInitial(value="<<value<<")");

    _useInitial = value;
} // useInitial


// ---------------------------------------------------------------------------------------------------------------------
// Get flag associated with using initial value term in time history expression.
bool
pylith::bc::NeumannTimeDependent::useInitial(void) const {
    return _useInitial;
} // useInitial


// ---------------------------------------------------------------------------------------------------------------------
// Use rate value term in time history expression.
void
pylith::bc::NeumannTimeDependent::useRate(const bool value) {
    PYLITH_COMPONENT_DEBUG("useRate(value="<<value<<")");

    _useRate = value;
} // useRate


// ---------------------------------------------------------------------------------------------------------------------
// Get flag associated with using rate value term in time history expression.
bool
pylith::bc::NeumannTimeDependent::useRate(void) const {
    return _useRate;
} // useRate


// ---------------------------------------------------------------------------------------------------------------------
// Use time history term in time history expression.
void
pylith::bc::NeumannTimeDependent::useTimeHistory(const bool value) {
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
    PYLITH_COMPONENT_DEBUG("setScaleName(value"<<value<<")");

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
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    pylith::feassemble::IntegratorBoundary* integrator = new pylith::feassemble::IntegratorBoundary(this);assert(integrator);
    integrator->setMarkerLabel(getMarkerLabel());
    integrator->setSubfieldName(getSubfieldName());
    integrator->setLabelName(getMarkerLabel());

    _NeumannTimeDependent::setKernelsResidual(integrator, *this, solution, _formulation);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
pylith::feassemble::Constraint*
pylith::bc::NeumannTimeDependent::createConstraint(const pylith::topology::Field& solution) {
    PYLITH_COMPONENT_DEBUG("createConstraint(solution="<<solution.getLabel()<<") empty method");

    return NULL;
} // createConstraint


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::bc::NeumannTimeDependent::createAuxiliaryField(const pylith::topology::Field& solution,
                                                       const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->setLabel("NeumannTimeDependent auxiliary field");

    assert(_auxiliaryFactory);
    assert(_normalizer);
    pylith::topology::Field::Description description = solution.getSubfieldInfo(_subfieldName.c_str()).description;
    if (_scaleName == std::string("pressure")) {
        description.scale = _normalizer->getPressureScale();
    } else if (_scaleName == std::string("velocity")) {
        description.scale = sqrt(_normalizer->getPressureScale() / _normalizer->getDensityScale());
    } else if (_scaleName == std::string("length")) {
        description.scale = _normalizer->getLengthScale();
    } else if (_scaleName == std::string("time")) {
        description.scale = sqrt(_normalizer->getDensityScale() / _normalizer->getPressureScale()) * _normalizer->getLengthScale();
    } else if (_scaleName == std::string("density")) {
        description.scale = _normalizer->getDensityScale();
    } else {
        std::ostringstream msg;
        msg << "Unknown name of scale ("<<_scaleName<<") for Neumann boundary condition for '" << _boundaryLabel << "'.";
        PYLITH_COMPONENT_ERROR(msg.str());
        throw std::logic_error(msg.str());
    } // if/else
    _auxiliaryFactory->initialize(auxiliaryField, *_normalizer, solution.getSpaceDim(), &description);

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
        if (_dbTimeHistory) {
            _dbTimeHistory->open();
        } // if
    } // _useTimeHistory

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    assert(_auxiliaryFactory);
    _auxiliaryFactory->setValuesFromDB();

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG("Displaying auxiliary field");
        auxiliaryField->view("Neumann auxiliary field");
    } // if

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::bc::NeumannTimeDependent::createDerivedField(const pylith::topology::Field& solution,
                                                     const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<") empty method");

    PYLITH_METHOD_RETURN(NULL);
} // createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary subfields at beginning of time step.
void
pylith::bc::NeumannTimeDependent::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                                       const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateAuxiliaryField(auxiliaryField="<<auxiliaryField<<", t="<<t<<")");

    if (_useTimeHistory) {
        assert(_normalizer);
        const PylithScalar timeScale = _normalizer->getTimeScale();
        TimeDependentAuxiliaryFactory::updateAuxiliaryField(auxiliaryField, t, timeScale, _dbTimeHistory);
    } // if

    PYLITH_METHOD_END;
} // updateAuxiliaryField


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
    PYLITH_COMPONENT_DEBUG("_updateKernelConstants(dt="<<dt<<")");

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
// Set kernels for residual.
void
pylith::bc::_NeumannTimeDependent::setKernelsResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                                      const pylith::bc::NeumannTimeDependent& bc,
                                                      const topology::Field& solution,
                                                      const pylith::problems::Physics::FormulationEnum formulation) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_NeumannTimeDependent::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "setKernelsResidual(integrator="<<integrator<<", bc="<<typeid(bc).name()<<", solution="
          << solution.getLabel()<<")"
          << pythia::journal::endl;

    const pylith::topology::Field::VectorFieldEnum fieldType = solution.getSubfieldInfo(bc.getSubfieldName()).description.vectorFieldType;
    const bool isScalarField = fieldType == pylith::topology::Field::SCALAR;

    const int bitInitial = bc.useInitial() ? 0x1 : 0x0;
    const int bitRate = bc.useRate() ? 0x2 : 0x0;
    const int bitTimeHistory = bc.useTimeHistory() ? 0x4 : 0x0;
    const int bitUse = bitInitial | bitRate | bitTimeHistory;

    PetscBdPointFunc r0 = NULL;
    PetscBdPointFunc r1 = NULL;
    switch (bitUse) {
    case 0x1:
        r0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::g0_initial_scalar :
             pylith::fekernels::NeumannTimeDependent::g0_initial_vector;
        break;
    case 0x2:
        r0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::g0_rate_scalar :
             pylith::fekernels::NeumannTimeDependent::g0_rate_vector;
        break;
    case 0x4:
        r0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::g0_timeHistory_scalar :
             pylith::fekernels::NeumannTimeDependent::g0_timeHistory_vector;
        break;
    case 0x3:
        r0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::g0_initialRate_scalar :
             pylith::fekernels::NeumannTimeDependent::g0_initialRate_vector;
        break;
    case 0x5:
        r0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::g0_initialTimeHistory_scalar :
             pylith::fekernels::NeumannTimeDependent::g0_initialTimeHistory_vector;
        break;
    case 0x6:
        r0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::g0_rateTimeHistory_scalar :
             pylith::fekernels::NeumannTimeDependent::g0_rateTimeHistory_vector;
        break;
    case 0x7:
        r0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::g0_initialRateTimeHistory_scalar :
             pylith::fekernels::NeumannTimeDependent::g0_initialRateTimeHistory_vector;
        break;
    case 0x0: {
        pythia::journal::warning_t warning(_NeumannTimeDependent::pyreComponent);
        warning << pythia::journal::at(__HERE__)
                << "Neumann time-dependent BC provides no values."
                << pythia::journal::endl;
        break;
    } // case 0x0
    default: {
        PYLITH_JOURNAL_LOGICERROR("Unknown combination of flags for Neumann BC terms (useInitial="
                                  <<bc.useInitial()
                                  << ", useRate="<<bc.useRate()<<", useTimeHistory="<<bc.useTimeHistory()<<").");
    } // default
    } // switch

    std::vector<ResidualKernels> kernels(1);
    switch (formulation) {
    case pylith::problems::Physics::QUASISTATIC:
        kernels[0] = ResidualKernels(bc.getSubfieldName(), pylith::feassemble::Integrator::RESIDUAL_LHS, r0, r1);
        break;
    case pylith::problems::Physics::DYNAMIC_IMEX:
    case pylith::problems::Physics::DYNAMIC:
        kernels[0] = ResidualKernels(bc.getSubfieldName(), pylith::feassemble::Integrator::RESIDUAL_RHS, r0, r1);
        break;
    default: {
        PYLITH_JOURNAL_LOGICERROR("Unknown formulation for equations ("<<formulation<<").");
    } // default
    } // switch

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // setKernelsResidual


// End of file
