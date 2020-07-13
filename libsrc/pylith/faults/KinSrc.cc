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
    _auxiliaryField(NULL),
    _slipLocalVec(NULL),
    _originTime(0.0) {} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrc::~KinSrc(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::KinSrc::deallocate(void) {
    PetscErrorCode err = VecDestroy(&_slipLocalVec);PYLITH_CHECK_ERROR(err);
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

    // Set default discretization of auxiliary subfields to match slip subfield in integrator auxiliary field.
    assert(_auxiliaryFactory);
    const pylith::topology::FieldBase::Discretization& discretization = faultAuxField.subfieldInfo("slip").fe;
    _auxiliaryFactory->setSubfieldDiscretization("default", discretization.basisOrder, discretization.quadOrder,
                                                 discretization.dimension, discretization.cellBasis, discretization.isBasisContinuous,
                                                 discretization.feSpace);

    delete _auxiliaryField;_auxiliaryField = new pylith::topology::Field(faultAuxField.mesh());assert(_auxiliaryField);
    _auxiliaryField->setLabel("kinsrc auxiliary");
    _auxiliaryFieldSetup(normalizer, cs);
    _auxiliaryField->subfieldsSetup();
    _auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(faultAuxField, *_auxiliaryField);
    _auxiliaryField->allocate();
    _auxiliaryField->zeroLocal();

    _auxiliaryFactory->setValuesFromDB();

    journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG("Displaying kinematic earthquake source auxiliary field");
        _auxiliaryField->view("KinSrc auxiliary field");
    } // if

    PetscErrorCode err = DMCreateLocalVector(faultAuxField.dmMesh(), &_slipLocalVec);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // initialize


// ----------------------------------------------------------------------
// Set slip subfield in fault integrator's auxiliary field at time t.
void
pylith::faults::KinSrc::updateSlip(pylith::topology::Field* const faultAuxField,
                                   const PylithScalar t,
                                   const PylithScalar timeScale) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateSlip(auxiliaryField="<<faultAuxField<<", t="<<t<<", timeScale="<<timeScale<<")");

    if (!_slipFnKernel || (t < _originTime)) {
        PYLITH_METHOD_END;
    } // if

    assert(faultAuxField);
    assert(_auxiliaryField);

    _setFEConstants(*faultAuxField); // Constants are attached to the auxiliary field.

    // Set update kernel for each fault auxiliary subfield.
    const pylith::string_vector& subfieldNames = faultAuxField->subfieldNames();
    const size_t numSubfields = subfieldNames.size();
    PetscPointFunc* subfieldKernels = (numSubfields > 0) ? new PetscPointFunc[numSubfields] : NULL;
    // By default, set all auxiliary subfield update kernels to NULL.
    for (size_t i = 0; i < numSubfields; ++i) {
        subfieldKernels[i] = NULL;
    } // for
    subfieldKernels[faultAuxField->subfieldInfo("slip").index] = _slipFnKernel;

    // Create local vector for slip for this source.
    PetscErrorCode err = 0;
    err = VecSet(_slipLocalVec, 0.0);PYLITH_CHECK_ERROR(err);

    PetscDM dmFaultAux = faultAuxField->dmMesh();
    err = PetscObjectCompose((PetscObject) dmFaultAux, "dmAux", (PetscObject) _auxiliaryField->dmMesh());PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmFaultAux, "A", (PetscObject) _auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);
    err = DMProjectFieldLocal(dmFaultAux, t, _slipLocalVec, subfieldKernels, INSERT_VALUES,
                              _slipLocalVec);PYLITH_CHECK_ERROR(err);
    delete[] subfieldKernels;subfieldKernels = NULL;

    // :TODO: @brad Replace this with operation on section array. Also remove call to zeroLocal() in
    // FaultCohesive::updateAuxiliaryField().

    // Add contribution of slip for this rupture to slip for all ruptures.
    err = VecAXPY(faultAuxField->localVector(), 1.0, _slipLocalVec);

    PYLITH_METHOD_END;
} // slip


// ----------------------------------------------------------------------
// Set constants used in finite-element integrations.
void
pylith::faults::KinSrc::_setFEConstants(const pylith::topology::Field& faultAuxField) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEConstants(faultAuxField="<<faultAuxField.getLabel()<<")");

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a NULL label) is
    // the correct one.
    PetscDS prob = NULL;
    PetscDM dmAux = faultAuxField.dmMesh();assert(dmAux);
    PetscErrorCode err = DMGetDS(dmAux, &prob);PYLITH_CHECK_ERROR(err);assert(prob);

    // Pointwise functions have been set in DS
    const int numConstants = 1;
    PylithScalar constants[numConstants];
    constants[0] = _originTime;
    err = PetscDSSetConstants(prob, numConstants, constants);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEConstants


// End of file
