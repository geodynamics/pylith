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
#include "pylith/feassemble/Integrator.hh" // USES Integrator

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <typeinfo> // USES typeid()

#include "pylith/materials/Material.hh" \
    // TEMPORARY

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputIntegrator::OutputIntegrator(pylith::feassemble::Integrator* const integrator) :
    _integrator(integrator) {
    PyreComponent::setName("outputintegrator");

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
// Verify configuration is acceptable.
void
pylith::meshio::OutputIntegrator::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::verifyConfiguration(solution="<<solution.label()<<")");

    assert(_integrator);

    const pylith::topology::Field* auxiliaryField = _integrator->getAuxiliaryField();
    const pylith::topology::Field* derivedField = _integrator->getDerivedField();

    // Info fields should be in integrator's auxiliary field.
    const size_t numInfoFields = _infoFields.size();
    if (!auxiliaryField && (numInfoFields > 0)) {
        std::ostringstream msg;
        msg << "Integrator has no auxiliary field, but output of information fields requested.\n"
            << "Information fields requested:";
        for (size_t i = 0; i < numInfoFields; ++i) {
            msg << "    " << _infoFields[i] << "\n";
        } // for
        throw std::runtime_error(msg.str());
    } else if ((numInfoFields > 0) && (std::string("all") != _infoFields[0])) {
        assert(auxiliaryField);
        for (size_t i = 0; i < numInfoFields; i++) {
            if (!auxiliaryField->hasSubfield(_infoFields[i].c_str())) {
                std::ostringstream msg;
                msg << "Could not find field '" << _infoFields[i] << "' in auxiliary field '" << auxiliaryField->label() << "' for output.";
                throw std::runtime_error(msg.str());
            } // if
        } // for
    } // if/else

    // Data fields should be in solution, integrator's auxiliary field, or integrator's derived field.
    const size_t numDataFields = _dataFields.size();
    if ((numDataFields > 0) && (std::string("all") != _dataFields[0])) {
        for (size_t i = 0; i < numDataFields; i++) {
            if (solution.hasSubfield(_dataFields[i].c_str())) { continue;}
            if (auxiliaryField && auxiliaryField->hasSubfield(_dataFields[i].c_str())) { continue;}
            if (derivedField && derivedField->hasSubfield(_dataFields[i].c_str())) { continue;}

            std::ostringstream msg;
            msg << "Could not find field '" << _dataFields[i] << "' in solution '" << solution.label()
                << "' or auxiliary field '" << (auxiliaryField ? auxiliaryField->label() : "NULL") << "' or derived field "
                << (derivedField ? derivedField->label() : "NULL") << "' for output.";
            throw std::runtime_error(msg.str());
        } // for
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Write output for step in solution.
void
pylith::meshio::OutputIntegrator::_writeInfo(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputIntegrator::_writeInfo()");

    if (!_integrator) { PYLITH_METHOD_END;}

    assert(_integrator);

    // :KLUDGE: Temporary code to set _label and _labelId if integrator ISA Material. This will go away when each
    // material has its own PetscDM.
    const pylith::materials::Material* const material = dynamic_cast<const pylith::materials::Material* const>(_integrator);
    if (material) {
        _temporarySetLabel("material-id", material->getMaterialId());
        PYLITH_COMPONENT_DEBUG("Setting OutputIntegrator label='material-id' and label id="<<material->getMaterialId()<<".");
    } // if

    const pylith::topology::Field* auxiliaryField = _integrator->getAuxiliaryField();
    if (!auxiliaryField) { PYLITH_METHOD_END;}

    const pylith::string_vector& infoNames = _infoNamesExpanded(auxiliaryField);

    const bool isInfo = true;
    const pylith::topology::Mesh& domainMesh = _integrator->getIntegrationDomainMesh();
    _open(domainMesh, isInfo);
    _openDataStep(0.0, domainMesh);

    const size_t numInfoFields = infoNames.size();
    for (size_t i = 0; i < numInfoFields; i++) {
        if (auxiliaryField->hasSubfield(infoNames[i].c_str())) {
            pylith::topology::Field* fieldBuffer = _getBuffer(*auxiliaryField, infoNames[i].c_str());assert(fieldBuffer);
            _appendField(0.0, fieldBuffer, domainMesh);
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
pylith::meshio::OutputIntegrator::_writeDataStep(const PylithReal t,
                                                 const PylithInt tindex,
                                                 const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputIntegrator::_writeDataStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.label()<<")");

    assert(_integrator);
    const pylith::topology::Field* auxiliaryField = _integrator->getAuxiliaryField();
    const pylith::topology::Field* derivedField = _integrator->getDerivedField();
    const pylith::topology::Mesh& domainMesh = _integrator->getIntegrationDomainMesh();

    const pylith::string_vector& dataNames = _dataNamesExpanded(solution, auxiliaryField, derivedField);

    _openDataStep(t, domainMesh);

    const size_t numDataFields = dataNames.size();
    for (size_t i = 0; i < numDataFields; i++) {
        if (solution.hasSubfield(dataNames[i].c_str())) {
            pylith::topology::Field* fieldBuffer = _getBuffer(solution, dataNames[i].c_str());assert(fieldBuffer);
            _appendField(t, fieldBuffer, domainMesh);
        } else if (auxiliaryField && auxiliaryField->hasSubfield(dataNames[i].c_str())) {
            pylith::topology::Field* fieldBuffer = _getBuffer(*auxiliaryField, dataNames[i].c_str());assert(fieldBuffer);
            _appendField(t, fieldBuffer, domainMesh);
        } else if (derivedField && derivedField->hasSubfield(dataNames[i].c_str())) {
            pylith::topology::Field* fieldBuffer = _getBuffer(*derivedField, dataNames[i].c_str());assert(fieldBuffer);
            _appendField(t, fieldBuffer, domainMesh);
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
