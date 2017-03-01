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

#include "OutputSolnNew.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <typeinfo> // USES typeid()

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnNew::OutputSolnNew(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSolnNew::~OutputSolnNew(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputSolnNew::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    OutputManagerNew::deallocate();

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set names of solution fields to output.
void
pylith::meshio::OutputSolnNew::vertexDataFields(const char* names[],
                                                const int numNames)
{ // vertexDataFields
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputSolnNew::vertexDataFields(names="<<names<<", numNames="<<numNames<<")");

    assert((names && numNames) || (!names && !numNames));

    _vertexDataFields.resize(numNames);
    for (int i=0; i < numNames; ++i) {
        assert(names[i]);
        _vertexDataFields[i] = names[i];
    } // for

    PYLITH_METHOD_END;
} // vertexDataFields

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::meshio::OutputSolnNew::verifyConfiguration(const pylith::topology::Field& solution) const
{ // verifyConfiguration
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputSolnNew::verifyConfiguration(solution="<<solution.label()<<")");

    const size_t numFields = _vertexDataFields.size();
    if (numFields > 0 && std::string("all") != _vertexDataFields[0]) {
        for (size_t iField=0; iField < numFields; iField++) {
            if (!solution.hasSubfield(_vertexDataFields[iField].c_str())) {
                std::ostringstream msg;
                msg << "Could not find field '" << _vertexDataFields[iField] << "' in solution '" << solution.label() << "' for output.";
                throw std::runtime_error(msg.str());
            } // if
        } // for
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration

// ----------------------------------------------------------------------
// Write solution at time step.
void
pylith::meshio::OutputSolnNew::writeTimeStep(const PylithReal t,
                                             const PylithInt timeStep,
                                             const pylith::topology::Field& solution)
{ // writeTimeStep
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputSolnNew::writeTimeStep(t="<<t<<", timeStep="<<timeStep<<", solution="<<solution.label()<<")");

    if (!this->shouldWrite(t, timeStep)) {
        PYLITH_METHOD_END;
    } // if

    const size_t numFields = _vertexDataFields.size();
    if (1 == numFields && std::string("all") == _vertexDataFields[0]) {
        PYLITH_JOURNAL_ERROR(":TODO: @brad Implement writing all subfields.");
    } else {
        for (size_t iField=0; iField < numFields; iField++) {
            if (!solution.hasSubfield(_vertexDataFields[iField].c_str())) {
                std::ostringstream msg;
                msg << "Could not find field '" << _vertexDataFields[iField] << "' in solution for output.";
                throw std::runtime_error(msg.str());
            } // if
            pylith::topology::Field& fieldBuffer = this->getBuffer(solution, _vertexDataFields[iField].c_str());
            this->appendVertexField(t, fieldBuffer, fieldBuffer.mesh());
        } // if/else
    } // if

    PYLITH_METHOD_END;
} // writeTimeStep


// End of file
