// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/faults/KinSrcAuxiliaryFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "pylith/scales/Scales.hh" // USES Scales

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*
#include "pylith/utils/Exceptions.hh" // USES Exception

#include <cassert>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::KinSrcAuxiliaryFactory::KinSrcAuxiliaryFactory(void) {
    GenericComponent::setName("kinsrcauxiliaryfactory");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrcAuxiliaryFactory::~KinSrcAuxiliaryFactory(void) {}


// ------------------------------------------------------------------------------------------------
// Add slip initiation time (relative to origin time) subfield to auxiliary field.
void
pylith::faults::KinSrcAuxiliaryFactory::addInitiationTime(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_DEBUG(pylith::journal::application_flow, "addInitiationTime(void)");

    const char* subfieldName = "initiation_time";
    const PylithReal timeScale = _scales->getTimeScale();

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


// ------------------------------------------------------------------------------------------------
// Add rise time subfield to auxiliary field.
void
pylith::faults::KinSrcAuxiliaryFactory::addRiseTime(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_DEBUG(pylith::journal::application_flow, "addRiseTime(void)");

    const char* subfieldName = "rise_time";
    const PylithReal timeScale = _scales->getTimeScale();

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


// ------------------------------------------------------------------------------------------------
// Add impulse duration subfield to auxiliary field.
void
pylith::faults::KinSrcAuxiliaryFactory::addImpulseDuration(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_DEBUG(pylith::journal::application_flow, "addImpulseDuration(void)");

    const char* subfieldName = "impulse_duration";
    const PylithReal timeScale = _scales->getTimeScale();

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
} // addImpulseDuration


// ------------------------------------------------------------------------------------------------
// Add final slip subfield to auxiliary field.
void
pylith::faults::KinSrcAuxiliaryFactory::addFinalSlip(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_DEBUG(pylith::journal::application_flow, "addFinalSlip(void)");

    const char* subfieldName = "final_slip";
    const char* componentNames[3] = { "final_slip_opening", "final_slip_left_lateral", "final_slip_reverse" };

    const PylithReal displacementScale = _scales->getDisplacementScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = displacementScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addFinalSlip


// ------------------------------------------------------------------------------------------------
// Add slip rate subfield to auxiliary field.
void
pylith::faults::KinSrcAuxiliaryFactory::addSlipRate(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_DEBUG(pylith::journal::application_flow, "addSlipRate(void)");

    const char* subfieldName = "slip_rate";
    const char* componentNames[3] = { "slip_rate_opening", "slip_rate_left_lateral", "slip_rate_reverse" };

    const PylithReal displacementScale = _scales->getDisplacementScale();
    const PylithReal timeScale = _scales->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = displacementScale / timeScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addSlipRate


// ------------------------------------------------------------------------------------------------
// Add time history value subfield to auxiliary field.
void
pylith::faults::KinSrcAuxiliaryFactory::addTimeHistoryValue(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_DEBUG(pylith::journal::application_flow, "addTimeHistoryValue(void)");

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
    // No subfield query; populated at beginning of time step.

    PYLITH_METHOD_END;
} // addTimeHistoryValue


// ------------------------------------------------------------------------------------------------
void
pylith::faults::KinSrcAuxiliaryFactory::updateTimeHistoryValue(pylith::topology::Field* auxiliaryField,
                                                               const PylithReal t,
                                                               const PylithReal timeScale,
                                                               spatialdata::spatialdb::TimeHistory* const dbTimeHistory) {
    PYLITH_METHOD_BEGIN;
    PYLITH_DEBUG(pylith::journal::application_flow_detail5, "Updating time history value, t="<<t<<", timeScale="<<timeScale);
    assert(auxiliaryField);
    assert(dbTimeHistory);

    PetscSection auxiliaryFieldSection = auxiliaryField->getLocalSection();assert(auxiliaryFieldSection);
    PetscInt pStart = 0, pEnd = 0;
    PylithCallPetsc(PetscSectionGetChart(auxiliaryFieldSection, &pStart, &pEnd));
    pylith::topology::VecVisitorMesh auxiliaryFieldVisitor(*auxiliaryField);
    PetscScalar* auxiliaryFieldArray = auxiliaryFieldVisitor.localArray();
    if (pEnd > pStart) { assert(auxiliaryFieldArray); }

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
                PYLITH_ERROR(pylith::ValueError, pylith::journal::external,
                             "Error querying for time '" << tDim << "' in time history database '" << dbTimeHistory->getDescription() << "'.");
            } // if
        } // if

        // Update value (normalized amplitude) in auxiliary field.
        const PetscInt offValue = auxiliaryFieldVisitor.sectionSubfieldOffset(i_value, p);
        auxiliaryFieldArray[offValue] = value;
    } // for

    PYLITH_METHOD_END;
} // updateAuilixaryField


// End of file
