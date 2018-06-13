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

#include "OutputMaterial.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/IntegratorPointwise.hh" // USES IntegratorPointwise

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <typeinfo> // USES typeid()

#include "pylith/materials/Material.hh" // TEMPORARY
#include <iostream>

// ----------------------------------------------------------------------
const char* pylith::meshio::OutputMaterial::_pyreComponent = "outputmaterial";

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputMaterial::OutputMaterial(pylith::feassemble::IntegratorPointwise* const integrator) :
    _integrator(integrator)
{ // constructor
    PyreComponent::name(_pyreComponent);

    // :KLUDGE: Temporary code to set _label and _labelId if integrator ISA Material.
    const pylith::materials::Material* const material = dynamic_cast<const pylith::materials::Material* const>(integrator);
    if (material) {
        _temporarySetLabel("material-id", material->id());
        PYLITH_COMPONENT_DEBUG("Setting OutputMaterial label='material-id' and label id="<<material->id()<<".");
    } // if
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputMaterial::~OutputMaterial(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputMaterial::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    OutputManager::deallocate();

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::meshio::OutputMaterial::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::verifyConfiguration(solution="<<solution.label()<<")");

    assert(_integrator);
    const pylith::topology::Field* auxField = _integrator->auxField();
    const pylith::topology::Field* derivedField = _integrator->derivedField();

    // Info fields should be in integrator's auxiliary field.
    const size_t numInfoFields = _infoFields.size();
    if (!auxField && (numInfoFields > 0)) {
        std::ostringstream msg;
        msg << "Integrator has no auxiliary field, but output of information fields requested.\n"
            << "Information fields requested:";
        for (size_t i = 0; i < numInfoFields; ++i) {
            msg << "    " << _infoFields[i] << "\n";
        } // for
        throw std::runtime_error(msg.str());
    } else if ((numInfoFields > 0) && (std::string("all") != _infoFields[0])) {
        assert(auxField);
        for (size_t i = 0; i < numInfoFields; i++) {
            if (!auxField->hasSubfield(_infoFields[i].c_str())) {
                std::ostringstream msg;
                msg << "Could not find field '" << _infoFields[i] << "' in auxiliary field '" << auxField->label() << "' for output.";
                throw std::runtime_error(msg.str());
            } // if
        } // for
    } // if/else

    // Data fields should be in solution, integrator's auxiliary field, or integrator's derived field.
    const size_t numDataFields = _dataFields.size();
    if ((numDataFields > 0) && (std::string("all") != _dataFields[0])) {
        for (size_t i = 0; i < numDataFields; i++) {
            if (solution.hasSubfield(_dataFields[i].c_str())) { continue; }
            if (auxField && auxField->hasSubfield(_dataFields[i].c_str())) { continue; }
            if (derivedField && derivedField->hasSubfield(_dataFields[i].c_str())) { continue; }

            std::ostringstream msg;
            msg << "Could not find field '" << _dataFields[i] << "' in solution '" << solution.label()
                << "' or auxiliary field '" << (auxField ? auxField->label() : "NULL") << "' or derived field "
                << (derivedField ? derivedField->label() : "NULL") << "' for output.";
            throw std::runtime_error(msg.str());
        } // for
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration

// ----------------------------------------------------------------------
// Write output for step in solution.
void
pylith::meshio::OutputMaterial::_writeInfo(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputMaterial::_writeInfo()");

    if (!_integrator) { PYLITH_METHOD_END; }

    assert(_integrator);
    const pylith::topology::Field* auxField = _integrator->auxField();
    if (!auxField) { PYLITH_METHOD_END; }

    const pylith::string_vector& infoNames = _infoNamesExpanded(auxField);

    const bool isInfo = true;
    _open(auxField->mesh(), isInfo);
    _openDataStep(0.0, auxField->mesh());

    const size_t numInfoFields = infoNames.size();
    for (size_t i = 0; i < numInfoFields; i++) {
        if (auxField->hasSubfield(infoNames[i].c_str())) {
            pylith::topology::Field* fieldBuffer = _getBuffer(*auxField, infoNames[i].c_str()); assert(fieldBuffer);
            _appendField(0.0, fieldBuffer, fieldBuffer->mesh());
        } else {
            std::ostringstream msg;
            msg << "Could not find field '" << infoNames[i] << "' in auxiliary field for info output.";
            throw std::runtime_error(msg.str());
        } // if/else
    } // for

    _closeDataStep();
    _close();

    PYLITH_METHOD_END;
} // _writeInfo


// ----------------------------------------------------------------------
// Write output for step in solution.
void
pylith::meshio::OutputMaterial::_writeDataStep(const PylithReal t,
                                               const PylithInt tindex,
                                               const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputMaterial::_writeDataStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.label()<<")");

    assert(_integrator);
    const pylith::topology::Field* auxField = _integrator->auxField();
    const pylith::topology::Field* derivedField = _integrator->derivedField();

    const pylith::string_vector& dataNames = _dataNamesExpanded(solution, auxField, derivedField);

    _openDataStep(t, solution.mesh());

    const size_t numDataFields = dataNames.size();
    for (size_t i = 0; i < numDataFields; i++) {
        if (solution.hasSubfield(dataNames[i].c_str())) {
            pylith::topology::Field* fieldBuffer = _getBuffer(solution, dataNames[i].c_str()); assert(fieldBuffer);
            _appendField(t, fieldBuffer, fieldBuffer->mesh());
        } else if (auxField && auxField->hasSubfield(dataNames[i].c_str())) {
            pylith::topology::Field* fieldBuffer = _getBuffer(*auxField, dataNames[i].c_str()); assert(fieldBuffer);
            _appendField(t, fieldBuffer, fieldBuffer->mesh());
        } else if (derivedField && derivedField->hasSubfield(dataNames[i].c_str())) {
            pylith::topology::Field* fieldBuffer = _getBuffer(*derivedField, dataNames[i].c_str()); assert(fieldBuffer);
            _appendField(t, fieldBuffer, fieldBuffer->mesh());
        } else {
            std::ostringstream msg;
            msg << "Could not find field '" << dataNames[i] << "' for data output.";
            throw std::runtime_error(msg.str());
        } // if/else
    } // for

    _closeDataStep();

    PYLITH_METHOD_END;
} // _writeDataStep


// End of file
