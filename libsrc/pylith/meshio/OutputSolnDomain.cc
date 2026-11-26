// -*- c++ -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


#include <portinfo>

#include "pylith/meshio/OutputSolnDomain.hh" // implementation of class methods

#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/exceptions/error.hh" // USES PYLITH_METHOD_*
#include "pylith/petsc/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

namespace pylith {
    namespace meshio {
        class _OutputSolnDomain {
          public:
            // Logging events
            class Events {
              public:
                static void init(void);

                static pylith::petsc::EventLogger logger;
                static PetscLogEvent writeSolnStep;
            }; // class Events
        }; // class _OutputSolnDomain
    } // namespace meshio
} // namespace pylith


// ------------------------------------------------------------------------------------------------
pylith::petsc::EventLogger pylith::meshio::_OutputSolnDomain::Events::logger;
PetscLogEvent pylith::meshio::_OutputSolnDomain::Events::writeSolnStep;


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_OutputSolnDomain::Events::init(void) {
    logger.setClassName("OutputSolnDomain");
    logger.initialize();
    writeSolnStep = logger.registerEvent("PL:OutputSolnDomain:writeSolnStep");
} // init


// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnDomain::OutputSolnDomain(void) {
    PyreComponent::setName("outputsolndomain");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSolnDomain::~OutputSolnDomain(void) {
}


// ---------------------------------------------------------------------------------------------------------------------
// Write data for step in solution.
void
pylith::meshio::OutputSolnDomain::_writeSolnStep(const PylithReal t,
                                                 const PylithInt tindex,
                                                 const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow,
                           "_writeSolnStep(t=" << t << ", tindex=" << tindex << ", solution=" << solution.getLabel()
                                               << ")");
    pylith::petsc::EventLogger::Scope logScope(_OutputSolnDomain::Events::writeSolnStep);

    const pylith::string_vector& subfieldNames = pylith::topology::FieldOps::getSubfieldNamesDomain(solution);
    PetscVec solutionVector = solution.getOutputVector();
    assert(solutionVector);

    const size_t numSubfieldNames = subfieldNames.size();
    for (size_t iField = 0; iField < numSubfieldNames; iField++) {
        assert(solution.hasSubfield(subfieldNames[iField].c_str()));

        OutputSubfield* subfield = NULL;
        subfield = OutputObserver::_getSubfield(solution, solution.getMesh(), subfieldNames[iField].c_str());
        assert(subfield);
        subfield->project(solutionVector);

        if (0 == iField) {
            // Need output mesh from subfield (which may be refined).
            assert(subfield);
            pylith::topology::Mesh* outputMesh = _getOutputMesh(*subfield);
            _openSolnStep(t, *outputMesh);
        } // if
        OutputObserver::_appendField(t, *subfield);
    } // for
    _closeSolnStep();

    PYLITH_METHOD_END;
} // _writeSolnStep


// End of file
