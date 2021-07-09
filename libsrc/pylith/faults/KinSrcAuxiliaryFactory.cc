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

#include "KinSrcAuxiliaryFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::KinSrcAuxiliaryFactory::KinSrcAuxiliaryFactory(void) {
    GenericComponent::setName("kinsrcauxiliaryfactory");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrcAuxiliaryFactory::~KinSrcAuxiliaryFactory(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add slip initiation time (relative to origin time) subfield to auxiliary field.
void
pylith::faults::KinSrcAuxiliaryFactory::addInitiationTime(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addInitiationTime(void)");

    const char* subfieldName = "initiation_time";
    const PylithReal timeScale = _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = timeScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addInitiationTime


// ---------------------------------------------------------------------------------------------------------------------
// Add riseTime subfield to auxiliary field.
void
pylith::faults::KinSrcAuxiliaryFactory::addRiseTime(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addRiseTime(void)");

    const char* subfieldName = "rise_time";
    const PylithReal timeScale = _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = timeScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addRiseTime


// ---------------------------------------------------------------------------------------------------------------------
// Add final slip subfield to auxiliary field.
void
pylith::faults::KinSrcAuxiliaryFactory::addFinalSlip(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addFinalSlip(void)");

    const char* subfieldName = "final_slip";
    const char* componentNames[3] = { "final_slip_opening", "final_slip_left_lateral", "final_slip_reverse" };

    const PylithReal lengthScale = _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = lengthScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addFinalSlip


// ---------------------------------------------------------------------------------------------------------------------
// Add slip rate subfield to auxiliary field.
void
pylith::faults::KinSrcAuxiliaryFactory::addSlipRate(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSlipRate(void)");

    const char* subfieldName = "slip_rate";
    const char* componentNames[3] = { "slip_rate_opening", "slip_rate_left_lateral", "slip_rate_reverse" };

    const PylithReal lengthScale = _normalizer->getLengthScale();
    const PylithReal timeScale = _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = lengthScale / timeScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addSlipRate


// ---------------------------------------------------------------------------------------------------------------------
// Add time history value subfield to auxiliary field.
void
pylith::faults::KinSrcAuxiliaryFactory::addTimeHistoryValue(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTimeHistoryValue(void)");

    const char* subfieldName = "time_history_value";

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization("final_slip"));
    // No subfield query; populated at begining of time step.

    PYLITH_METHOD_END;
} // addTimeHistoryValue


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::faults::KinSrcAuxiliaryFactory::updateTimeHistoryValue(pylith::topology::Field* auxiliaryField,
                                                               const PylithReal t,
                                                               const PylithReal timeScale,
                                                               spatialdata::spatialdb::TimeHistory* const dbTimeHistory) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug("kinsrcauxiliaryfactory");
    debug << pythia::journal::at(__HERE__)
          << "KinSrcAuxiliaryFactory::updateTimeHistoryValue(auxiliaryField="<<auxiliaryField<<", t="<<t
          <<", timeScale="<<timeScale<<", dbTimeHistory="<<dbTimeHistory<<")"
          << pythia::journal::endl;

    assert(auxiliaryField);
    assert(dbTimeHistory);

    PetscErrorCode err = 0;

    PetscSection auxiliaryFieldSection = auxiliaryField->getLocalSection();assert(auxiliaryFieldSection);
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryFieldSection, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    pylith::topology::VecVisitorMesh auxiliaryFieldVisitor(*auxiliaryField);
    PetscScalar* auxiliaryFieldArray = auxiliaryFieldVisitor.localArray();assert(auxiliaryFieldArray);

    // Compute offset of time history subfields in auxiliary field.
    const PetscInt i_startTime = auxiliaryField->getSubfieldInfo("initiation_time").index;
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
                msg << "Error querying for time '" << tDim << "' in time history database '" << dbTimeHistory->getLabel() << "'.";
                throw std::runtime_error(msg.str());
            } // if
        } // if

        // Update value (normalized amplitude) in auxiliary field.
        const PetscInt offValue = auxiliaryFieldVisitor.sectionSubfieldOffset(i_value, p);
        auxiliaryFieldArray[offValue] = value;
    } // for

    PYLITH_METHOD_END;
} // updateAuilixaryField


// End of file
