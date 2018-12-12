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

#include "OutputPhysics.hh" // Implementation of class methods

#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/meshio/FieldFilter.hh" // USES FieldFilter
#include "pylith/meshio/OutputTrigger.hh" // USES OutputTrigger
#include "pylith/feassemble/PhysicsImplementation.hh" // USES PhysicsImplementation

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <iostream> // USES std::cout
#include <typeinfo> // USES typeid()

#include "pylith/materials/Material.hh" // TEMPORARY

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputPhysics::OutputPhysics(void) :
    _label(""),
    _labelId(0) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputPhysics::~OutputPhysics(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputPhysics::deallocate(void) {
    ObserverPhysics::deallocate();
    OutputObserver::deallocate();
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set names of vertex information fields to output.
void
pylith::meshio::OutputPhysics::setInfoFields(const char* names[],
                                             const int numNames) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::setInfoFields(names="<<names<<", numNames="<<numNames<<")");

    assert((names && numNames) || (!names && !numNames));

    _infoFieldNames.resize(numNames);
    for (int i = 0; i < numNames; ++i) {
        assert(names[i]);
        _infoFieldNames[i] = names[i];
    } // for

    PYLITH_METHOD_END;
} // setInfoFields


// ---------------------------------------------------------------------------------------------------------------------
// Get names of vertex information fields to output.
const pylith::string_vector&
pylith::meshio::OutputPhysics::getInfoFields(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::setInfoFields()");

    PYLITH_METHOD_RETURN(_infoFieldNames);
} // getInfoFields


// ---------------------------------------------------------------------------------------------------------------------
// Set names of vertex data fields to output.
void
pylith::meshio::OutputPhysics::setDataFields(const char* names[],
                                             const int numNames) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::setDataFields(names="<<names<<", numNames="<<numNames<<")");

    assert((names && numNames) || (!names && !numNames));

    _dataFieldNames.resize(numNames);
    for (int i = 0; i < numNames; ++i) {
        assert(names[i]);
        _dataFieldNames[i] = names[i];
    } // for

    PYLITH_METHOD_END;
} // setDataFields


// ---------------------------------------------------------------------------------------------------------------------
// Get names of vertex data fields to output.
const pylith::string_vector&
pylith::meshio::OutputPhysics::getDataFields(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::getDataFields()");

    PYLITH_METHOD_RETURN(_dataFieldNames);
} // getDataFields


// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::meshio::OutputPhysics::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::verifyConfiguration(solution="<<solution.label()<<")");

    assert(_physics);

    const pylith::topology::Field* auxiliaryField = _physics->getAuxiliaryField();
    const pylith::topology::Field* derivedField = _physics->getDerivedField();

    // Info fields should be in auxiliary field.
    const size_t numInfoFields = _infoFieldNames.size();
    if (!auxiliaryField && (numInfoFields > 0)) {
        std::ostringstream msg;
        msg << "Physics implementation '" << _physics->getName() << "' has no auxiliary field, but physics output '"
            << PyreComponent::getIdentifier() << "' requested information fields:\n";
        for (size_t i = 0; i < numInfoFields; ++i) {
            msg << "    " << _infoFieldNames[i] << "\n";
        } // for
        throw std::runtime_error(msg.str());
    } else if ((numInfoFields > 0) && (std::string("all") != _infoFieldNames[0])) {
        assert(auxiliaryField);
        for (size_t i = 0; i < numInfoFields; i++) {
            if (!auxiliaryField->hasSubfield(_infoFieldNames[i].c_str())) {
                std::ostringstream msg;
                msg << "Could not find subfield '" << _infoFieldNames[i] << "' in auxiliary field '"
                    << auxiliaryField->label() << "' for physics output '" << PyreComponent::getIdentifier() << "''.";
                throw std::runtime_error(msg.str());
            } // if
        } // for
    } // if/else

    // Data fields should be in solution, constraint's auxiliary field, or constraint's derived field.
    const size_t numDataFields = _dataFieldNames.size();
    if ((numDataFields > 0) && (std::string("all") != _dataFieldNames[0])) {
        for (size_t i = 0; i < numDataFields; i++) {
            if (solution.hasSubfield(_dataFieldNames[i].c_str())) { continue;}
            if (auxiliaryField && auxiliaryField->hasSubfield(_dataFieldNames[i].c_str())) { continue;}
            if (derivedField && derivedField->hasSubfield(_dataFieldNames[i].c_str())) { continue;}

            std::ostringstream msg;
            msg << "Could not find subfield '" << _dataFieldNames[i] << "' in solution field '" << solution.label()
                << ", auxiliary field '" << (auxiliaryField ? auxiliaryField->label() : "NULL") << "', or derived field "
                << (derivedField ? derivedField->label() : "NULL") << "' for physics output '"
                << PyreComponent::getIdentifier() << "'.";
            throw std::runtime_error(msg.str());
        } // for
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Get update from integrator (subject of observer).
void
pylith::meshio::OutputPhysics::update(const PylithReal t,
                                      const PylithInt tindex,
                                      const pylith::topology::Field& solution,
                                      const bool infoOnly) {
    if (infoOnly) {
        _writeInfo();
    } else {
        assert(_trigger);
        if (_trigger->shouldWrite(t, tindex)) {
            _writeDataStep(t, tindex, solution);
        } // if
    } // if/else
} // update


// ----------------------------------------------------------------------
// Write output for step in solution.
void
pylith::meshio::OutputPhysics::_writeInfo(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::_writeInfo()");

    if (!_physics) { PYLITH_METHOD_END;}

    assert(_physics);

    // :KLUDGE: Temporary code to set _label and _labelId if physics ISA Material. This will go away when each
    // material has its own PetscDM.
    const pylith::materials::Material* const material = dynamic_cast<const pylith::materials::Material* const>(_physics);
    if (material) {
        _temporarySetLabel("material-id", material->getMaterialId());
        PYLITH_COMPONENT_DEBUG("Setting label='material-id' and label id="<<material->getMaterialId()<<".");
    } // if

    const pylith::topology::Field* auxiliaryField = _physics->getAuxiliaryField();
    if (!auxiliaryField) { PYLITH_METHOD_END;}

    const pylith::string_vector& infoNames = _expandInfoFieldNames(auxiliaryField);

    const bool isInfo = true;
    const pylith::topology::Mesh& domainMesh = _physics->getPhysicsDomainMesh();
    _open(domainMesh, isInfo);
    _openDataStep(0.0, domainMesh);

    const size_t numInfoFields = infoNames.size();
    for (size_t i = 0; i < numInfoFields; i++) {
        if (auxiliaryField->hasSubfield(infoNames[i].c_str())) {
            pylith::topology::Field* fieldBuffer = _getBuffer(*auxiliaryField, infoNames[i].c_str());assert(fieldBuffer);
            _appendField(0.0, fieldBuffer, domainMesh);
        } else {
            std::ostringstream msg;
            msg << "Internal Error: Could not find subfield '" << infoNames[i] << "' in auxiliary field for info output.";
            PYLITH_COMPONENT_ERROR(msg.str());
            throw std::runtime_error(msg.str());
        } // if/else
    } // for

    _closeDataStep();
    _close();

    PYLITH_METHOD_END;
} // _writeInfo


// ---------------------------------------------------------------------------------------------------------------------
// Prepare for output.
void
pylith::meshio::OutputPhysics::_open(const pylith::topology::Mesh& mesh,
                                     const bool isInfo) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::open(mesh="<<typeid(mesh).name()<<", isInfo="<<isInfo<<")");

    if (!_writer) {
        PYLITH_COMPONENT_ERROR("Writer for physics output '" << PyreComponent::getIdentifier() << "' not set.");
    } // if

    assert(_writer);
    _writer->open(mesh, isInfo, _label.length() ? _label.c_str() : NULL, _labelId);

    PYLITH_METHOD_END;
} // _open


// ---------------------------------------------------------------------------------------------------------------------
/// Close output files.
void
pylith::meshio::OutputPhysics::_close(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::_close()");

    assert(_writer);
    _writer->close();

    PYLITH_METHOD_END;
} // _close


// ---------------------------------------------------------------------------------------------------------------------
// Prepare for output at this solution step.
void
pylith::meshio::OutputPhysics::_openDataStep(const PylithReal t,
                                             const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::_openDataStep(t="<<t<<", mesh="<<typeid(mesh).name()<<")");

    assert(_writer);
    if (!_writer->isOpen()) {
        bool infoOnly = false;
        _writer->open(mesh, infoOnly, _label.length() ? _label.c_str() : NULL, _labelId);
    } // if
    _writer->openTimeStep(t, mesh, _label.length() ? _label.c_str() : NULL, _labelId);

    PYLITH_METHOD_END;
} // _openDataStep


// ---------------------------------------------------------------------------------------------------------------------
// Finalize output at this solution step.
void
pylith::meshio::OutputPhysics::_closeDataStep(void) {
    assert(_writer);
    _writer->closeTimeStep();

} // _closeDataStep


// ----------------------------------------------------------------------
// Write output for step in solution.
void
pylith::meshio::OutputPhysics::_writeDataStep(const PylithReal t,
                                              const PylithInt tindex,
                                              const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::_writeDataStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.label()<<")");

    assert(_physics);
    const pylith::topology::Field* auxiliaryField = _physics->getAuxiliaryField();
    const pylith::topology::Field* derivedField = _physics->getDerivedField();
    const pylith::topology::Mesh& domainMesh = _physics->getPhysicsDomainMesh();

    const pylith::string_vector& dataNames = _expandDataFieldNames(solution, auxiliaryField, derivedField);

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
            msg << "Internal Error: Could not find subfield '" << dataNames[i] << "' for data output.";
            PYLITH_COMPONENT_ERROR(msg.str());
            throw std::runtime_error(msg.str());
        } // if/else
    } // for

    _closeDataStep();

    PYLITH_METHOD_END;
} // _writeDataStep


// ---------------------------------------------------------------------------------------------------------------------
// Append finite-element vertex field to file.
void
pylith::meshio::OutputPhysics::_appendField(const PylithReal t,
                                            pylith::topology::Field* field,
                                            const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::appendField(t="<<t<<", field="<<typeid(field).name()<<", mesh="<<typeid(mesh).name()<<")");

    assert(field);

    pylith::topology::Field* fieldFiltered = _fieldFilter->filter(field);
    pylith::topology::Field* fieldDimensioned = _dimensionField(fieldFiltered);assert(fieldDimensioned);

    const int basisOrder = _getBasisOrder(*fieldDimensioned);
    switch (basisOrder) {
    case 0:
        _writer->writeCellField(t, *fieldDimensioned, _label.c_str(), _labelId);
        break;

    case 1:
        _writer->writeVertexField(t, *fieldDimensioned, mesh);
        break;

    default:
        PYLITH_COMPONENT_ERROR(
            "Unsupported basis order for output ("
                << basisOrder <<"). Use FieldFilterProject with basis order of 0 or 1. Skipping output of '"
                << field->label() << "' field."
            );
    } // switch

    PYLITH_METHOD_END;
} // _appendField


// ---------------------------------------------------------------------------------------------------------------------
// Names of information fields for output.
pylith::string_vector
pylith::meshio::OutputPhysics::_expandInfoFieldNames(const pylith::topology::Field* auxField) const {
    PYLITH_METHOD_BEGIN;

    if (auxField && (1 == _infoFieldNames.size()) && (std::string("all") == _infoFieldNames[0])) {
        PYLITH_METHOD_RETURN(auxField->subfieldNames());
    } // if

    PYLITH_METHOD_RETURN(_infoFieldNames);
} // _expandInfoFieldNames


// ---------------------------------------------------------------------------------------------------------------------
// Names of data fields for output.
pylith::string_vector
pylith::meshio::OutputPhysics::_expandDataFieldNames(const pylith::topology::Field& solution,
                                                     const pylith::topology::Field* auxField,
                                                     const pylith::topology::Field* derivedField) const {
    PYLITH_METHOD_BEGIN;

    if ((1 == _dataFieldNames.size()) && (std::string("all") == _dataFieldNames[0])) {
        pylith::string_vector dataNames;
        dataNames = solution.subfieldNames();

        if (derivedField) {
            const pylith::string_vector& derivedSubfields = derivedField->subfieldNames();
            const size_t numAdd = derivedSubfields.size();
            dataNames.resize(dataNames.size() + numAdd);
            for (size_t iAdd = 0, iName = dataNames.size(); iAdd < numAdd; ++iAdd) {
                dataNames[iName] = derivedSubfields[iAdd];
            } // for
        } // if
        PYLITH_METHOD_RETURN(dataNames);
    } // if

    PYLITH_METHOD_RETURN(_dataFieldNames);
} // _expandDataFieldNames


// ---------------------------------------------------------------------------------------------------------------------
// TEMPOARY Set label and label id.
void
pylith::meshio::OutputPhysics::_temporarySetLabel(const char* label,
                                                  const PylithInt labelId) {
    _label = label;
    _labelId = labelId;
} // _temporarySetLabel


// End of file
