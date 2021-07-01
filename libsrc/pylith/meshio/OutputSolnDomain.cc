// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "OutputSolnDomain.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnDomain::OutputSolnDomain(void) {
    PyreComponent::setName("outputsolndomain");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSolnDomain::~OutputSolnDomain(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Write data for step in solution.
void
pylith::meshio::OutputSolnDomain::_writeSolnStep(const PylithReal t,
                                                 const PylithInt tindex,
                                                 const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_writeSolnStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.getLabel()<<")");

    const pylith::string_vector& subfieldNames = pylith::topology::FieldOps::getSubfieldNamesDomain(solution);
    PetscVec solutionVector = solution.getOutputVector();assert(solutionVector);

    _openSolnStep(t, solution.getMesh());
    const size_t numSubfieldNames = subfieldNames.size();
    for (size_t iField = 0; iField < numSubfieldNames; iField++) {
        assert(solution.hasSubfield(subfieldNames[iField].c_str()));

        OutputSubfield* subfield = NULL;
        subfield = OutputObserver::_getSubfield(solution, solution.getMesh(), subfieldNames[iField].c_str());assert(subfield);
        subfield->project(solutionVector);

        OutputObserver::_appendField(t, *subfield);
    } // for
    _closeSolnStep();

    PYLITH_METHOD_END;
} // _writeSolnStep


// End of file
