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

#include "OutputIntegrator.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/journals.hh" \
    // USES PYLITH_COMPONENT_*


// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputIntegrator::OutputIntegrator(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputIntegrator::~OutputIntegrator(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputIntegrator::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    OutputManager::deallocate();

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set names of information fields to output.
void
pylith::meshio::OutputIntegrator::vertexInfoFields(const char* names[],
                                                   const int numNames) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputIntegrator::vertexInfoFields(names="<<names<<", numNames="<<numNames<<")");

    assert((names && numNames) || (!names && !numNames));

    _vertexInfoFields.resize(numNames);
    for (int i = 0; i < numNames; ++i) {
        assert(names[i]);
        _vertexInfoFields[i] = names[i];
    } // for

    PYLITH_METHOD_END;
} // vertexInfoFields

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::meshio::OutputIntegrator::verifyConfiguration(const pylith::topology::Field& auxField) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputIntegrator::verifyConfiguration(auxField="<<auxField.label()<<")");

    const size_t numFields = _vertexInfoFields.size();
    if ((numFields > 0) && (std::string("all") != _vertexInfoFields[0])) {
        for (size_t iField = 0; iField < numFields; iField++) {
            if (!auxField.hasSubfield(_vertexInfoFields[iField].c_str())) {
                std::ostringstream msg;
                msg << "Could not find field '" << _vertexInfoFields[iField] << "' in auxiliary field '" << auxField.label() << "' for output.";
                throw std::runtime_error(msg.str());
            } // if
        } // for
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration

// ----------------------------------------------------------------------
// Write information.
void
pylith::meshio::OutputIntegrator::writeInfo(const pylith::topology::Field& auxField) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputIntegrator::writeTimeStep(auxField="<<auxField.label()<<")");

    const pylith::string_vector& subfieldNames = (1 == _vertexInfoFields.size() && std::string("all") == _vertexInfoFields[0]) ? auxField.subfieldNames() : _vertexInfoFields;

    const bool isInfo = true;
    this->open(auxField.mesh(), isInfo);
    this->openTimeStep(0.0, auxField.mesh());
    const size_t numFields = subfieldNames.size();
    for (size_t iField = 0; iField < numFields; iField++) {
        if (!auxField.hasSubfield(subfieldNames[iField].c_str())) {
            std::ostringstream msg;
            msg << "Could not find field '" << subfieldNames[iField] << "' in auxiliary field for output.";
            throw std::runtime_error(msg.str());
        } // if

        pylith::topology::Field& fieldBuffer = this->getBuffer(auxField, subfieldNames[iField].c_str());
        this->appendVertexField(0.0, fieldBuffer, fieldBuffer.mesh());
    } // for
    this->closeTimeStep();
    this->close();

    PYLITH_METHOD_END;
} // writeInfo


// ----------------------------------------------------------------------
// Write solution at time step.
void
pylith::meshio::OutputIntegrator::writeTimeStep(const PylithReal t,
                                                const PylithInt tindex,
                                                const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputMaterial::writeTimeStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.label()<<") empty method");
    PYLITH_METHOD_END;
} // writeTimeStep


// End of file
