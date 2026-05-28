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

#include "pylith/bc/TimeDependentOps.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/scales/Scales.hh" // USES Nondimensional

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/Exceptions.hh" // USES Exceptions

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
void
pylith::bc::TimeDependentOps::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                                   const pylith::real t,
                                                   const pylith::real timeScale,
                                                   const std::shared_ptr<spatialdata::spatialdb::TimeHistory>& dbTimeHistory) {
    PYLITH_METHOD_BEGIN;
    PYLITH_DEBUG(pylith::journal::application_flow_detail5, "Updating auxiliary field.");
    assert(auxiliaryField);
    assert(dbTimeHistory);

    PetscSection auxiliaryFieldSection = auxiliaryField->getLocalSection();assert(auxiliaryFieldSection);
    pylith::integer pStart = 0, pEnd = 0;
    PylithCallPetsc(PetscSectionGetChart(auxiliaryFieldSection, &pStart, &pEnd));
    if (pStart == pEnd) {
        PYLITH_METHOD_END;
    } // if

    pylith::topology::VecVisitorMesh auxiliaryFieldVisitor(*auxiliaryField);
    pylith::scalar* auxiliaryFieldArray = auxiliaryFieldVisitor.localArray();assert(auxiliaryFieldArray);

    // Compute offset of time history subfields in auxiliary field.
    const pylith::integer i_startTime = auxiliaryField->getSubfieldInfo("time_history_start_time").index;
    const pylith::integer i_value = auxiliaryField->getSubfieldInfo("time_history_value").index;

    // Loop over all points in section.
    for (pylith::integer p = pStart; p < pEnd; ++p) {
        // Skip points without values in section.
        if (!auxiliaryFieldVisitor.sectionDof(p)) {continue;}

        // Get starting time and compute relative time for point.
        const pylith::integer offStartTime = auxiliaryFieldVisitor.sectionSubfieldOffset(i_startTime, p);
        const pylith::real tStart = auxiliaryFieldArray[offStartTime];
        const pylith::real tRel = t - tStart;

        // Query time history for value (normalized amplitude).
        pylith::scalar value = 0.0;
        if (tRel >= 0.0) {
            pylith::real tDim = tRel * timeScale;
            const int err = dbTimeHistory->query(&value, tDim);
            if (err) {
                PYLITH_COMPONENT_ERROR(pylith::ValueError, pylith::journal::external,
                                       "Error querying for time '" << tDim << "' in time history database '" << dbTimeHistory->getDescription() << "'.");
            } // if
        } // if

        // Update value (normalized amplitude) in auxiliary field.
        const pylith::integer offValue = auxiliaryFieldVisitor.sectionSubfieldOffset(i_value, p);
        auxiliaryFieldArray[offValue] = value;
    } // for

    PYLITH_METHOD_END;
} // updateAuxiliaryField
