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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "OutputSoln.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <typeinfo> // USES typeid()

// ----------------------------------------------------------------------
const char* pylith::meshio::OutputSoln::_pyreComponent = "outputsoln";

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSoln::OutputSoln(void) {
    PyreComponent::name(_pyreComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSoln::~OutputSoln(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputSoln::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    OutputManager::deallocate();

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::meshio::OutputSoln::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputSoln::verifyConfiguration(solution="<<solution.label()<<")");

    const size_t numSubfields = _dataFields.size();
    if ((numSubfields > 0) && (std::string("all") != _dataFields[0])) {
        for (size_t iField = 0; iField < numSubfields; iField++) {
            if (!solution.hasSubfield(_dataFields[iField].c_str())) {
                std::ostringstream msg;
                msg << "Could not find field '" << _dataFields[iField] << "' in solution '" << solution.label() << "' for output.";
                throw std::runtime_error(msg.str());
            } // if
        } // for
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration

// ----------------------------------------------------------------------
// Write data for step in solution.
void
pylith::meshio::OutputSoln::_writeDataStep(const PylithReal t,
                                           const PylithInt tindex,
                                           const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputSoln::_writeDataStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.label()<<")");

    const pylith::string_vector& subfieldNames = (1 == _dataFields.size() && std::string("all") == _dataFields[0]) ? solution.subfieldNames() : _dataFields;

    _openDataStep(t, solution.mesh());
    const size_t numSubfields = subfieldNames.size();
    for (size_t iField = 0; iField < numSubfields; iField++) {
        if (!solution.hasSubfield(subfieldNames[iField].c_str())) {
            std::ostringstream msg;
            msg << "Could not find field '" << subfieldNames[iField] << "' in solution for output.";
            throw std::runtime_error(msg.str());
        } // if

        pylith::topology::Field* fieldBuffer = _getBuffer(solution, subfieldNames[iField].c_str()); assert(fieldBuffer);
        _appendField(t, *fieldBuffer, fieldBuffer->mesh());
    } // for
    _closeDataStep();

    PYLITH_METHOD_END;
} // _writeDataStep


// End of file
