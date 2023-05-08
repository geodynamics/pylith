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

#include "TractionPerturbation.hh" // implementation of object methods

#include "pylith/bc/TimeDependentAuxiliaryFactory.hh" // USES TimeDependentAuxiliaryFactory

#include "pylith/fekernels/TimeDependentFn.hh" // USES TimeDependentFn kernels

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace faults {
        class _TractionPerturbation {
            // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////
public:

            static const char* pyreComponent;

        }; // _TractionPerturbation
        const char* _TractionPerturbation::pyreComponent = "tractionperturbation";

    } // bc
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::TractionPerturbation::TractionPerturbation(void) :
    _dbTimeHistory(NULL),
    _auxiliaryFactory(new pylith::bc::TimeDependentAuxiliaryFactory(pylith::bc::TimeDependentAuxiliaryFactory::TANGENTIAL_NORMAL)),
    _auxiliaryField(NULL),
    _tractionKernel(NULL),
    _useInitial(true),
    _useRate(false),
    _useTimeHistory(false),
    _originTime(0.0) {
    PyreComponent::setName(_TractionPerturbation::pyreComponent);
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::TractionPerturbation::~TractionPerturbation(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::TractionPerturbation::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
    _dbTimeHistory = NULL; // :KLUDGE: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set time history database.
void
pylith::faults::TractionPerturbation::setTimeHistoryDB(spatialdata::spatialdb::TimeHistory* th) {
    PYLITH_COMPONENT_DEBUG("setTimeHistoryDB(th"<<th<<")");

    _dbTimeHistory = th;
} // setTimeHistoryDB


// ------------------------------------------------------------------------------------------------
// Get time history database.
const spatialdata::spatialdb::TimeHistory*
pylith::faults::TractionPerturbation::getTimeHistoryDB(void) {
    return _dbTimeHistory;
} // getTimeHistoryDB


// ------------------------------------------------------------------------------------------------
// Use initial value term in time history expression.
void
pylith::faults::TractionPerturbation::useInitial(const bool value) {
    PYLITH_COMPONENT_DEBUG("useInitial(value="<<value<<")");

    _useInitial = value;
} // useInitial


// ------------------------------------------------------------------------------------------------
// Get flag associated with using initial value term in time history expression.
bool
pylith::faults::TractionPerturbation::useInitial(void) const {
    return _useInitial;
} // useInitial


// ------------------------------------------------------------------------------------------------
// Use rate value term in time history expression.
void
pylith::faults::TractionPerturbation::useRate(const bool value) {
    PYLITH_COMPONENT_DEBUG("useRate(value="<<value<<")");

    _useRate = value;
} // useRate


// ------------------------------------------------------------------------------------------------
// Get flag associated with using rate value term in time history expression.
bool
pylith::faults::TractionPerturbation::useRate(void) const {
    return _useRate;
} // useRate


// ------------------------------------------------------------------------------------------------
// Use time history term in time history expression.
void
pylith::faults::TractionPerturbation::useTimeHistory(const bool value) {
    PYLITH_COMPONENT_DEBUG("useTimeHistory(value="<<value<<")");

    _useTimeHistory = value;
} // useTimeHistory


// ------------------------------------------------------------------------------------------------
// Get flag associated with using time history term in time history expression.
bool
pylith::faults::TractionPerturbation::useTimeHistory(void) const { // useTimeHistory
    return _useTimeHistory;
} // useTimeHistory


// ------------------------------------------------------------------------------------------------
// Set origin time for earthquake source.
void
pylith::faults::TractionPerturbation::setOriginTime(const PylithReal value) {
    _originTime = value;
}


// ------------------------------------------------------------------------------------------------
// Get origin time for earthquake source.
PylithReal
pylith::faults::TractionPerturbation::getOriginTime(void) const {
    return _originTime;
}


// ------------------------------------------------------------------------------------------------
// Get auxiliary field associated with the kinematic source.
const pylith::topology::Field&
pylith::faults::TractionPerturbation::getAuxiliaryField(void) const {
    PYLITH_METHOD_BEGIN;

    assert(_auxiliaryField);
    PYLITH_METHOD_RETURN(*_auxiliaryField);
} // getAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Set the spatial database for filling auxiliary subfields.
void
pylith::faults::TractionPerturbation::setAuxiliaryFieldDB(spatialdata::spatialdb::SpatialDB* value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setAuxiliaryFieldDB(value="<<typeid(value).name()<<")");

    assert(_auxiliaryFactory);
    _auxiliaryFactory->setQueryDB(value);

    PYLITH_METHOD_END;
} // setAuxiliaryFieldDB


// ------------------------------------------------------------------------------------------------
// Initialize kinematic (prescribed slip) earthquake source.
void
pylith::faults::TractionPerturbation::initialize(const pylith::topology::Field& faultAuxField,
                                                 const spatialdata::units::Nondimensional& normalizer,
                                                 const spatialdata::geocoords::CoordSys* cs) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(faultAuxField"<<faultAuxField.getLabel()<<", normalizer, cs="<<typeid(cs).name()<<")");

    // Set default discretization of auxiliary subfields to match slip/slip_rate subfield in integrator auxiliary field.
    assert(_auxiliaryFactory);
    const pylith::topology::FieldBase::Discretization& discretization = faultAuxField.getSubfieldInfo("tractionPerturbation").fe;
    _auxiliaryFactory->setSubfieldDiscretization("default", discretization.basisOrder, discretization.quadOrder,
                                                 discretization.dimension, discretization.isFaultOnly,
                                                 discretization.cellBasis, discretization.feSpace,
                                                 discretization.isBasisContinuous);

    delete _auxiliaryField;_auxiliaryField = new pylith::topology::Field(faultAuxField.getMesh());assert(_auxiliaryField);
    _auxiliaryField->setLabel("perturbation auxiliary field");
    _auxiliaryFieldSetup(normalizer, cs);
    pylith::topology::FieldOps::checkDiscretization(faultAuxField, *_auxiliaryField);
    _auxiliaryField->allocate();

    _auxiliaryFactory->setValuesFromDB();

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG("Displaying kinematic earthquake source auxiliary field");
        _auxiliaryField->view("KinSrc auxiliary field");
    } // if

    PYLITH_METHOD_END;
} // initialize


// ------------------------------------------------------------------------------------------------
// Get requested slip subfields at time t.
void
pylith::faults::TractionPerturbation::getTraction(PetscVec tractionLocalVec,
                                                  pylith::topology::Field* faultAuxiliaryField,
                                                  const PylithScalar t,
                                                  const PylithScalar timeScale) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getTraction(tractionLocalVec="<<tractionLocalVec<<", faultAuxiliaryField="<<faultAuxiliaryField
                                                          <<", t="<<t<<", timeScale="<<timeScale<<")");

    if (!_tractionKernel || (t < _originTime)) {
        PYLITH_METHOD_END;
    } // if

    assert(tractionLocalVec);
    assert(_auxiliaryField);

    // Update auxiliary field
    if (_useTimeHistory) {
        pylith::bc::TimeDependentAuxiliaryFactory::updateAuxiliaryField(_auxiliaryField, t, timeScale, _dbTimeHistory);
    } // if
    _setFEConstants(*faultAuxiliaryField); // Constants are attached to the auxiliary field for the traction vector.

    const size_t numSubfields = 2;
    assert(faultAuxiliaryField->getSubfieldNames().size() == numSubfields);
    PetscPointFunc subfieldKernels[numSubfields] = {
        _tractionKernel,
        NULL, // fault rheology traction
    };

    // Create local vector for traction for this perturbation.
    PetscErrorCode err = 0;
    PetscDM faultAuxiliaryDM = faultAuxiliaryField->getDM();
    PetscDMLabel dmLabel = NULL;
    PetscInt labelValue = 0;
    const PetscInt part = 0;
    err = DMSetAuxiliaryVec(faultAuxiliaryDM, dmLabel, labelValue, part,
                            _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMProjectFieldLocal(faultAuxiliaryDM, t-_originTime, tractionLocalVec, subfieldKernels, INSERT_VALUES,
                              tractionLocalVec);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // getTraction


// ------------------------------------------------------------------------------------------------
// Set kernels for computing traction perturbation.
void
pylith::faults::TractionPerturbation::_setTractionKernel(void) {
    PYLITH_METHOD_BEGIN;

    const int bitInitial = useInitial() ? 0x1 : 0x0;
    const int bitRate = useRate() ? 0x2 : 0x0;
    const int bitTimeHistory = useTimeHistory() ? 0x4 : 0x0;
    const int bitUse = bitInitial | bitRate | bitTimeHistory;

    switch (bitUse) {
    case 0x1:
        _tractionKernel = pylith::fekernels::TimeDependentFn::initial_vector;
        break;
    case 0x2:
        _tractionKernel = pylith::fekernels::TimeDependentFn::rate_vector;
        break;
    case 0x4:
        _tractionKernel = pylith::fekernels::TimeDependentFn::timeHistory_vector;
        break;
    case 0x3:
        _tractionKernel = pylith::fekernels::TimeDependentFn::initialRate_vector;
        break;
    case 0x5:
        _tractionKernel = pylith::fekernels::TimeDependentFn::initialTimeHistory_vector;
        break;
    case 0x6:
        _tractionKernel = pylith::fekernels::TimeDependentFn::rateTimeHistory_vector;
        break;
    case 0x7:
        _tractionKernel = pylith::fekernels::TimeDependentFn::initialRateTimeHistory_vector;
        break;
    case 0x0: {
        pythia::journal::warning_t warning(_TractionPerturbation::pyreComponent);
        warning << pythia::journal::at(__HERE__)
                << "Neumann time-dependent BC provides no values."
                << pythia::journal::endl;
        break;
    } // case 0x0
    default: {
        PYLITH_COMPONENT_LOGICERROR("Unknown combination of flags for Neumann BC terms (useInitial="
                                    <<useInitial() << ", useRate="<<useRate()<<", useTimeHistory="
                                    <<useTimeHistory()<<").");
    } // default
    } // switch

    PYLITH_METHOD_END;
} // setTractionKernel


// ------------------------------------------------------------------------------------------------
// Preinitialize earthquake source. Set names/sizes of auxiliary subfields.
void
pylith::faults::TractionPerturbation::_auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                                                           const spatialdata::geocoords::CoordSys* cs) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxiliaryFieldSetup()");

    assert(_auxiliaryFactory);
    assert(cs);
    _auxiliaryFactory->initialize(_auxiliaryField, normalizer, cs->getSpaceDim());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the slip time function
    // kernel.

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

    _auxiliaryField->subfieldsSetup();
    _auxiliaryField->createDiscretization();

    PYLITH_METHOD_END;
} // _auxiliaryFieldSetup


// End of file
