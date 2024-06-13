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

#include "AuxiliaryFactorySourceTime.hh" // implementation of object methods

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
pylith::sources::AuxiliaryFactorySourceTime::AuxiliaryFactorySourceTime(void) {
    GenericComponent::setName("auxiliaryfactorysourcetime");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::AuxiliaryFactorySourceTime::~AuxiliaryFactorySourceTime(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add center frequency of source time function to auxiliary fields.
void
pylith::sources::AuxiliaryFactorySourceTime::addCenterFrequency(void) { // centerFrequency
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
pylith::sources::AuxiliaryFactorySourceTime::addTimeHistoryStartTime(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTimeHistoryStartTime(void)");

    const char* subfieldName = "time_history_start_time";

    assert(_normalizer);
    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = _normalizer->getTimeScale();
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addTimeHistoryStartTime


// ---------------------------------------------------------------------------------------------------------------------
// Add time history value field to auxiliary fields.
void
pylith::sources::AuxiliaryFactorySourceTime::addTimeHistoryValue(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTimeHistoryValue(void)");

    const char* subfieldName = "time_history_value";

    assert(_normalizer);
    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization("time_history_amplitude"));
    // No subfield query; populated by integrator or constraint at beginning of time step.

    PYLITH_METHOD_END;
} // addTimeHistoryValue


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::sources::AuxiliaryFactorySourceTime::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
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
