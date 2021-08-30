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

#include "KinSrc.hh" // implementation of object methods

#include "pylith/faults/KinSrcAuxiliaryFactory.hh" // USES KinSrcAuxiliaryFactory
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps::checkDisretization()

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <typeinfo> // USES typeid()
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::KinSrc::KinSrc(void) :
    _auxiliaryFactory(new pylith::faults::KinSrcAuxiliaryFactory),
    _slipFnKernel(NULL),
    _slipRateFnKernel(NULL),
    _slipAccFnKernel(NULL),
    _auxiliaryField(NULL),
    _originTime(0.0) {}


// ----------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrc::~KinSrc(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::KinSrc::deallocate(void) {
    delete _auxiliaryField;_auxiliaryField = NULL;
    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ----------------------------------------------------------------------
// Set origin time for earthquake source.
void
pylith::faults::KinSrc::originTime(const PylithReal value) {
    _originTime = value;
} // originTime


// ----------------------------------------------------------------------
// Get origin time for earthquake source.
PylithReal
pylith::faults::KinSrc::originTime(void) const {
    return _originTime;
} // originTime


// ----------------------------------------------------------------------
// Get auxiliary field.
const pylith::topology::Field&
pylith::faults::KinSrc::auxField(void) const {
    PYLITH_METHOD_BEGIN;

    assert(_auxiliaryField);

    PYLITH_METHOD_RETURN(*_auxiliaryField);
} // auxField


// ----------------------------------------------------------------------
// Set database for auxiliary fields.
void
pylith::faults::KinSrc::auxFieldDB(spatialdata::spatialdb::SpatialDB* value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("auxFieldDB(value="<<typeid(value).name()<<")");

    assert(_auxiliaryFactory);
    _auxiliaryFactory->setQueryDB(value);

    PYLITH_METHOD_END;
} // auxFieldDB


// ----------------------------------------------------------------------
// Initialize kinematic (prescribed slip) earthquake source.
void
pylith::faults::KinSrc::initialize(const pylith::topology::Field& faultAuxField,
                                   const spatialdata::units::Nondimensional& normalizer,
                                   const spatialdata::geocoords::CoordSys* cs) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(faultAuxField"<<faultAuxField.getLabel()<<", normalizer, cs="<<typeid(cs).name()<<")");

    // Set default discretization of auxiliary subfields to match slip/slip_rate subfield in integrator auxiliary field.
    assert(_auxiliaryFactory);
    const char* slipFieldName = faultAuxField.hasSubfield("slip") ? "slip" : "slip_acceleration";
    const pylith::topology::FieldBase::Discretization& discretization = faultAuxField.getSubfieldInfo(slipFieldName).fe;
    _auxiliaryFactory->setSubfieldDiscretization("default", discretization.basisOrder, discretization.quadOrder,
                                                 discretization.dimension, discretization.isFaultOnly,
                                                 discretization.cellBasis, discretization.feSpace,
                                                 discretization.isBasisContinuous);

    delete _auxiliaryField;_auxiliaryField = new pylith::topology::Field(faultAuxField.getMesh());assert(_auxiliaryField);
    _auxiliaryField->setLabel("kinsrc auxiliary");
    _auxiliaryFieldSetup(normalizer, cs);
    _auxiliaryField->subfieldsSetup();
    _auxiliaryField->createDiscretization();
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


// ----------------------------------------------------------------------
// Set slip values at time t.
void
pylith::faults::KinSrc::updateSlip(PetscVec slipLocalVec,
                                   pylith::topology::Field* faultAuxiliaryField,
                                   const PylithScalar t,
                                   const PylithScalar timeScale) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateSlip(slipLocalVec="<<slipLocalVec<<", faultAuxiliaryField="<<faultAuxiliaryField
                                                     <<", t="<<t<<", timeScale="<<timeScale<<")");

    if (!_slipFnKernel || (t < _originTime)) {
        PYLITH_METHOD_END;
    } // if

    assert(slipLocalVec);
    assert(_auxiliaryField);

    _setFEConstants(*faultAuxiliaryField); // Constants are attached to the auxiliary field for the slip vector.

    PetscPointFunc subfieldKernels[1];
    subfieldKernels[0] = _slipFnKernel;

    // Create local vector for slip for this source.
    PetscErrorCode err = 0;
    PetscDM faultAuxiliaryDM = faultAuxiliaryField->getDM();
    PetscDMLabel dmLabel = NULL;
    PetscInt labelValue = 0;
    err = DMSetAuxiliaryVec(faultAuxiliaryDM, dmLabel, labelValue,
                            _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMProjectFieldLocal(faultAuxiliaryDM, t, slipLocalVec, subfieldKernels, INSERT_VALUES,
                              slipLocalVec);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // updateSlip


// ----------------------------------------------------------------------
// Set slip rate values at time t.
void
pylith::faults::KinSrc::updateSlipRate(PetscVec slipRateLocalVec,
                                       pylith::topology::Field* faultAuxiliaryField,
                                       const PylithScalar t,
                                       const PylithScalar timeScale) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateSlipRate(slipRateLocalVec="<<slipRateLocalVec<<", faultAuxiliaryField="<<faultAuxiliaryField
                                                             <<", t="<<t<<", timeScale="<<timeScale<<")");

    if (!_slipRateFnKernel || (t < _originTime)) {
        PYLITH_METHOD_END;
    } // if

    assert(slipRateLocalVec);
    assert(_auxiliaryField);

    _setFEConstants(*faultAuxiliaryField); // Constants are attached to the auxiliary field for the slip rate vector.

    PetscPointFunc subfieldKernels[1];
    subfieldKernels[0] = _slipRateFnKernel;

    // Create local vector for slip for this source.
    PetscErrorCode err = 0;
    PetscDM faultAuxiliaryDM = faultAuxiliaryField->getDM();
    PetscDMLabel dmLabel = NULL;
    PetscInt labelValue = 0;
    err = DMSetAuxiliaryVec(faultAuxiliaryDM, dmLabel, labelValue,
                            _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMProjectFieldLocal(faultAuxiliaryDM, t, slipRateLocalVec, subfieldKernels, INSERT_VALUES,
                              slipRateLocalVec);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // updateSlipRate


// ----------------------------------------------------------------------
// Set slip acceleration values at time t.
void
pylith::faults::KinSrc::updateSlipAcc(PetscVec slipAccLocalVec,
                                      pylith::topology::Field* faultAuxiliaryField,
                                      const PylithScalar t,
                                      const PylithScalar timeScale) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateSlipAcc(slipAccLocalVec="<<slipAccLocalVec<<", faultAuxiliaryField="<<faultAuxiliaryField
                                                           <<", t="<<t<<", timeScale="<<timeScale<<")");

    if (!_slipAccFnKernel || (t < _originTime)) {
        PYLITH_METHOD_END;
    } // if

    assert(slipAccLocalVec);
    assert(_auxiliaryField);

    _setFEConstants(*faultAuxiliaryField); // Constants are attached to the auxiliary field for the slip rate vector.

    PetscPointFunc subfieldKernels[1];
    subfieldKernels[0] = _slipAccFnKernel;

    // Create local vector for slip for this source.
    PetscErrorCode err = 0;
    PetscDM faultAuxiliaryDM = faultAuxiliaryField->getDM();
    PetscDMLabel dmLabel = NULL;
    PetscInt labelValue = 0;
    err = DMSetAuxiliaryVec(faultAuxiliaryDM, dmLabel, labelValue,
                            _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMProjectFieldLocal(faultAuxiliaryDM, t, slipAccLocalVec, subfieldKernels, INSERT_VALUES,
                              slipAccLocalVec);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // updateSlipAcc


// ----------------------------------------------------------------------
// Set constants used in finite-element integrations.
void
pylith::faults::KinSrc::_setFEConstants(const pylith::topology::Field& faultAuxField) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEConstants(faultAuxField="<<faultAuxField.getLabel()<<")");

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a NULL label) is
    // the correct one.
    PetscDS prob = NULL;
    PetscDM dmAux = faultAuxField.getDM();assert(dmAux);
    PetscErrorCode err = DMGetDS(dmAux, &prob);PYLITH_CHECK_ERROR(err);assert(prob);

    // Pointwise functions have been set in DS
    const int numConstants = 1;
    PylithScalar constants[numConstants];
    constants[0] = _originTime;
    err = PetscDSSetConstants(prob, numConstants, constants);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEConstants


// End of file
