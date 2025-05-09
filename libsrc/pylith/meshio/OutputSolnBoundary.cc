// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/meshio/OutputSolnBoundary.hh" // implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES createLowerDimMesh()
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnBoundary::OutputSolnBoundary(void) :
    _boundaryMesh(NULL),
    _labelName(""),
    _labelValue(1) {
    PyreComponent::setName("outputsolnboundary");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSolnBoundary::~OutputSolnBoundary(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputSolnBoundary::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    OutputSoln::deallocate();

    delete _boundaryMesh;_boundaryMesh = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set name of label identifier for subdomain.
void
pylith::meshio::OutputSolnBoundary::setLabelName(const char* value) {
    PYLITH_METHOD_BEGIN;

    _labelName = value;

    PYLITH_METHOD_END;
} // setLabelName


// ------------------------------------------------------------------------------------------------
// Set value of label identifier for subdomain.
void
pylith::meshio::OutputSolnBoundary::setLabelValue(const int value) {
    PYLITH_METHOD_BEGIN;

    _labelValue = value;

    PYLITH_METHOD_END;
} // setLabelValue


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::meshio::OutputSolnBoundary::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    OutputSoln::verifyConfiguration(solution);

    PetscDM dmSoln = solution.getDM();assert(dmSoln);
    PetscBool hasLabel = PETSC_FALSE;
    PetscErrorCode err = DMHasLabel(dmSoln, _labelName.c_str(), &hasLabel);PYLITH_CHECK_ERROR(err);
    if (!hasLabel) {
        std::ostringstream msg;
        msg << "Mesh missing group of points '" << _labelName << " for output using solution boundary observer '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
// Write data for step in solution.
void
pylith::meshio::OutputSolnBoundary::_writeSolnStep(const PylithReal t,
                                                   const PylithInt tindex,
                                                   const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_writeSolnStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.getLabel()<<")");

    if (!_boundaryMesh) {
        const char* componentName = this->getFullIdentifier();
        _boundaryMesh = pylith::topology::MeshOps::createLowerDimMesh(solution.getMesh(), _labelName.c_str(), _labelValue, componentName);
        assert(_boundaryMesh);
    } // if

    const pylith::string_vector& subfieldNames = pylith::topology::FieldOps::getSubfieldNamesDomain(solution);
    PetscVec solutionVector = solution.getOutputVector();assert(solutionVector);

    const size_t numSubfieldNames = subfieldNames.size();
    for (size_t iField = 0; iField < numSubfieldNames; iField++) {
        assert(solution.hasSubfield(subfieldNames[iField].c_str()));

        OutputSubfield* subfield = NULL;
        subfield = OutputObserver::_getSubfield(solution, *_boundaryMesh, subfieldNames[iField].c_str());assert(subfield);
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
