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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "AuxiliaryFactoryMomentTensorForce.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::sources::AuxiliaryFactoryMomentTensorForce::AuxiliaryFactoryMomentTensorForce(void) {
    GenericComponent::setName("auxiliaryfactorymomenttensorforce");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::AuxiliaryFactoryMomentTensorForce::~AuxiliaryFactoryMomentTensorForce(void) {}


// ----------------------------------------------------------------------------
// Add isotropic permeability subfield to auxiliary fields.
void
pylith::sources::AuxiliaryFactoryMomentTensorForce::addMomentTensor(void) { // momentTensor
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addMomentTensor(void)");

    const char* subfieldName = "moment_tensor";
    const char* componentNames[6] = {
        "moment_tensor_xx",
        "moment_tensor_yy",
        "moment_tensor_zz",
        "moment_tensor_xy",
        "moment_tensor_yz",
        "moment_tensor_xz"
    };
    const int tensorSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = (3 == _spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = tensorSize;
    description.componentNames.resize(tensorSize);
    for (int i = 0; i < tensorSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addMomentTensor


// ---------------------------------------------------------------------------------------------------------------------
// Add time delay of source time function to auxiliary fields.
void
pylith::sources::AuxiliaryFactoryMomentTensorForce::addTimeDelay(void) { // timeDelay
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTimeDelay(void)");

    const char* subfieldName = "time_delay";

    pylith::topology::Field::Description description;
    const PylithReal timeScale = _normalizer->getTimeScale();
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = timeScale;
    // description.validator = pylith::topology::FieldQuery::validatorNonnegative;
    description.validator = NULL;


    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addTimeDelay


// ---------------------------------------------------------------------------------------------------------------------
// Add center frequency of source time function to auxiliary fields.
void
pylith::sources::AuxiliaryFactoryMomentTensorForce::addCenterFrequency(void) { // centerFrequency
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addCenterFrequency(void)");

    const char* subfieldName = "center_frequency";
    const PylithReal timeScale = _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = 1.0 / timeScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addCenterFrequency


// ---------------------------------------------------------------------------------------------------------------------
// Add time history start time field to auxiliary fields.
void
pylith::sources::AuxiliaryFactoryMomentTensorForce::addTimeHistoryStartTime(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTimeHistoryStartTime(void)");

    const char* subfieldName = "time_history_start_time";

    assert(_defaultDescription);
    assert(_normalizer);
    pylith::topology::FieldBase::Description subfieldDescription(*_defaultDescription);
    subfieldDescription.label = subfieldName;
    subfieldDescription.alias = subfieldName;
    subfieldDescription.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    subfieldDescription.numComponents = 1;
    subfieldDescription.componentNames.resize(1);
    subfieldDescription.componentNames[0] = subfieldName;
    subfieldDescription.scale = _normalizer->getTimeScale();
    subfieldDescription.validator = NULL;

    _field->subfieldAdd(subfieldDescription, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addTimeHistoryStartTime


// ---------------------------------------------------------------------------------------------------------------------
// Add time history value field to auxiliary fields.
void
pylith::sources::AuxiliaryFactoryMomentTensorForce::addTimeHistoryValue(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTimeHistoryValue(void)");

    const char* subfieldName = "time_history_value";

    assert(_defaultDescription);
    assert(_normalizer);
    pylith::topology::FieldBase::Description subfieldDescription(*_defaultDescription);
    subfieldDescription.label = subfieldName;
    subfieldDescription.alias = subfieldName;
    subfieldDescription.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    subfieldDescription.numComponents = 1;
    subfieldDescription.componentNames.resize(1);
    subfieldDescription.componentNames[0] = subfieldName;
    subfieldDescription.validator = NULL;

    _field->subfieldAdd(subfieldDescription, getSubfieldDiscretization("time_history_amplitude"));
    // No subfield query; populated by integrator or constraint at beginning of time step.

    PYLITH_METHOD_END;
} // addTimeHistoryValue


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::sources::AuxiliaryFactoryMomentTensorForce::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                                                         const PylithReal t,
                                                                         const PylithReal timeScale,
                                                                         spatialdata::spatialdb::TimeHistory* const dbTimeHistory) {
    PYLITH_METHOD_BEGIN;
    // pythia::journal::debug_t debug(_TimeDependentAuxiliaryFactory::genericComponent);
    // debug << pythia::journal::at(__HERE__)
    //       << "TimeDependentAuxiliaryFactory::updateAuxiliaryField(auxiliaryField="<<auxiliaryField<<", t="<<t
    //       <<", timeScale="<<timeScale<<", dbTimeHistory="<<dbTimeHistory<<")"
    //       << pythia::journal::endl;

    assert(auxiliaryField);
    assert(dbTimeHistory);

    PetscErrorCode err = 0;

    PetscSection auxiliaryFieldSection = auxiliaryField->getLocalSection();assert(auxiliaryFieldSection);
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryFieldSection, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    pylith::topology::VecVisitorMesh auxiliaryFieldVisitor(*auxiliaryField);
    PetscScalar* auxiliaryFieldArray = auxiliaryFieldVisitor.localArray();assert(auxiliaryFieldArray);

    // Compute offset of time history subfields in auxiliary field.
    const PetscInt i_startTime = auxiliaryField->getSubfieldInfo("time_history_start_time").index;
    const PetscInt i_value = auxiliaryField->getSubfieldInfo("time_history_value").index;

    // Loop over all points in section.
    for (PetscInt p = pStart; p < pEnd; ++p) {
        // Skip points without values in section.
        if (!auxiliaryFieldVisitor.sectionDof(p)) {continue;}

        // Get starting time and compute relative time for point.
        const PetscInt offStartTime = auxiliaryFieldVisitor.sectionSubfieldOffset(i_startTime, p);
        const PylithScalar tStart = auxiliaryFieldArray[offStartTime];
        const PylithScalar tRel = t - tStart;

        // Query time history for value (normalized amplitude).
        PylithScalar value = 0.0;
        if (tRel >= 0.0) {
            PylithScalar tDim = tRel * timeScale;
            const int err = dbTimeHistory->query(&value, tDim);
            if (err) {
                std::ostringstream msg;
                msg << "Error querying for time '" << tDim << "' in time history database '" << dbTimeHistory->getDescription() << "'.";
                throw std::runtime_error(msg.str());
            } // if
        } // if

        // Update value (normalized amplitude) in auxiliary field.
        const PetscInt offValue = auxiliaryFieldVisitor.sectionSubfieldOffset(i_value, p);
        auxiliaryFieldArray[offValue] = value;
    } // for

    PYLITH_METHOD_END;
} // updateAuxiliaryField


// End of file
