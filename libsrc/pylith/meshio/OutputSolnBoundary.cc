// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "OutputSolnBoundary.hh" // implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END


// ----------------------------------------------------------------------
const char* pylith::meshio::OutputSolnBoundary::_pyreComponent = "outputsoln";

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnBoundary::OutputSolnBoundary(pylith::problems::Problem* const problem) :
    OutputSoln(problem),
    _label(""),
    _boundaryMesh(NULL)
{ // constructor
    PyreComponent::name(_pyreComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSolnBoundary::~OutputSolnBoundary(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputSolnBoundary::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    OutputSoln::deallocate();

    delete _boundaryMesh; _boundaryMesh = NULL;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set label identifier for subdomain.
void
pylith::meshio::OutputSolnBoundary::label(const char* value) {
    PYLITH_METHOD_BEGIN;

    _label = value;

    PYLITH_METHOD_END;
} // label

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::meshio::OutputSolnBoundary::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    OutputSoln::verifyConfiguration(solution);

    PetscDM dmMesh = solution.dmMesh();assert(dmMesh);
    PetscBool hasLabel = PETSC_FALSE;
    PetscErrorCode err = DMHasLabel(dmMesh, _label.c_str(), &hasLabel);PYLITH_CHECK_ERROR(err);
    if (!hasLabel) {
        std::ostringstream msg;
        msg << "Mesh missing group of vertices '" << _label << " for subdomain output.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Write data for step in solution.
void
pylith::meshio::OutputSolnBoundary::_writeDataStep(const PylithReal t,
                                                   const PylithInt tindex,
                                                   const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_writeDataStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.label()<<")");

    if (!_boundaryMesh) {
        _boundaryMesh = new pylith::topology::Mesh(solution.mesh(), _label.c_str()); assert(_boundaryMesh);
    } // if

    const pylith::topology::Field* auxField = NULL;
    const pylith::topology::Field* derivedField = NULL;

    const pylith::string_vector& dataNames = _dataNamesExpanded(solution, auxField, derivedField);

    _openDataStep(t, *_boundaryMesh);
    const size_t numDataFields = dataNames.size();
    for (size_t iField = 0; iField < numDataFields; iField++) {
        if (!solution.hasSubfield(dataNames[iField].c_str())) {
            std::ostringstream msg;
            msg << "Could not find field '" << dataNames[iField] << "' in solution for output.";
            throw std::runtime_error(msg.str());
        } // if

        pylith::topology::Field* fieldBuffer = _getBuffer(solution, dataNames[iField].c_str()); assert(fieldBuffer);
        _appendField(t, fieldBuffer, *_boundaryMesh);
    } // for
    _closeDataStep();

    PYLITH_METHOD_END;
} // _writeDataStep


// End of file
