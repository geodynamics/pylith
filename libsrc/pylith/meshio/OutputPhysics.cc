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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "OutputPhysics.hh" // Implementation of class methods

#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/meshio/OutputTrigger.hh" // USES OutputTrigger
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield
#include "pylith/feassemble/PhysicsImplementation.hh" // USES PhysicsImplementation

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <iostream> // USES std::cout
#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputPhysics::OutputPhysics(void) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputPhysics::~OutputPhysics(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputPhysics::deallocate(void) {
    ObserverPhysics::deallocate();
    OutputObserver::deallocate();
} // deallocate


// ------------------------------------------------------------------------------------------------
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


// ------------------------------------------------------------------------------------------------
// Get names of vertex information fields to output.
const pylith::string_vector&
pylith::meshio::OutputPhysics::getInfoFields(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::setInfoFields()");

    PYLITH_METHOD_RETURN(_infoFieldNames);
} // getInfoFields


// ------------------------------------------------------------------------------------------------
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


// ------------------------------------------------------------------------------------------------
// Get names of vertex data fields to output.
const pylith::string_vector&
pylith::meshio::OutputPhysics::getDataFields(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::getDataFields()");

    PYLITH_METHOD_RETURN(_dataFieldNames);
} // getDataFields


// ------------------------------------------------------------------------------------------------
// Set time scale.
void
pylith::meshio::OutputPhysics::setTimeScale(const PylithReal value) {
    OutputObserver::setTimeScale(value);
} // setTimeScale


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::meshio::OutputPhysics::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::verifyConfiguration(solution="<<solution.getLabel()<<")");

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
                    << auxiliaryField->getLabel() << "' for physics output '" << PyreComponent::getIdentifier() << "''.";
                throw std::runtime_error(msg.str());
            } // if
        } // for
    } // if/else

    // Data fields should be in solution, auxiliary field, or derived field.
    const size_t numDataFields = _dataFieldNames.size();
    if ((numDataFields > 0) && (std::string("all") != _dataFieldNames[0])) {
        for (size_t i = 0; i < numDataFields; i++) {
            if (solution.hasSubfield(_dataFieldNames[i].c_str())) { continue;}
            if (auxiliaryField && auxiliaryField->hasSubfield(_dataFieldNames[i].c_str())) { continue;}
            if (derivedField && derivedField->hasSubfield(_dataFieldNames[i].c_str())) { continue;}

            std::ostringstream msg;
            msg << "Could not find subfield '" << _dataFieldNames[i] << "' in solution field '" << solution.getLabel()
                << ", auxiliary field '" << (auxiliaryField ? auxiliaryField->getLabel() : "NULL") << "', or derived field "
                << (derivedField ? derivedField->getLabel() : "NULL") << "' for physics output '"
                << PyreComponent::getIdentifier() << "'.";
            throw std::runtime_error(msg.str());
        } // for
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
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


// ------------------------------------------------------------------------------------------------
// Write output for step in solution.
void
pylith::meshio::OutputPhysics::_writeInfo(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::_writeInfo()");

    if (!_physics) { PYLITH_METHOD_END;}

    assert(_physics);
    const pylith::topology::Field* auxiliaryField = _physics->getAuxiliaryField();
    if (!auxiliaryField) { PYLITH_METHOD_END;}

    const pylith::string_vector& infoNames = _expandInfoFieldNames(auxiliaryField);

    const bool isInfo = true;
    const pylith::topology::Mesh& domainMesh = _physics->getPhysicsDomainMesh();
    _open(domainMesh, isInfo);
    _openDataStep(0.0, domainMesh);

    PetscVec auxiliaryVector = auxiliaryField->getOutputVector();
    auxiliaryField->scatterLocalToOutput();

    const size_t numInfoFields = infoNames.size();
    for (size_t i = 0; i < numInfoFields; i++) {
        if (auxiliaryField->hasSubfield(infoNames[i].c_str())) {
            OutputSubfield* subfield = _getSubfield(*auxiliaryField, domainMesh, infoNames[i].c_str());
            subfield->project(auxiliaryVector);
            _appendField(0.0, *subfield);
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


// ------------------------------------------------------------------------------------------------
// Prepare for output.
void
pylith::meshio::OutputPhysics::_open(const pylith::topology::Mesh& mesh,
                                     const bool isInfo) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::open(mesh="<<typeid(mesh).name()<<", isInfo="<<isInfo<<")");

    if (!_writer) {
        PYLITH_COMPONENT_ERROR("Writer for physics output '" << PyreComponent::getIdentifier() << "' not set.");
    } // if

    assert(_trigger);
    _trigger->setTimeScale(_timeScale);

    assert(_writer);
    _writer->setTimeScale(_timeScale);
    _writer->open(mesh, isInfo);

    PYLITH_METHOD_END;
} // _open


// ------------------------------------------------------------------------------------------------
/// Close output files.
void
pylith::meshio::OutputPhysics::_close(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::_close()");

    assert(_writer);
    _writer->close();

    PYLITH_METHOD_END;
} // _close


// ------------------------------------------------------------------------------------------------
// Prepare for output at this solution step.
void
pylith::meshio::OutputPhysics::_openDataStep(const PylithReal t,
                                             const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::_openDataStep(t="<<t<<", mesh="<<typeid(mesh).name()<<")");

    assert(_writer);
    if (!_writer->isOpen()) {
        bool infoOnly = false;
        _open(mesh, infoOnly);
    } // if
    _writer->openTimeStep(t, mesh);

    PYLITH_METHOD_END;
} // _openDataStep


// ------------------------------------------------------------------------------------------------
// Finalize output at this solution step.
void
pylith::meshio::OutputPhysics::_closeDataStep(void) {
    assert(_writer);
    _writer->closeTimeStep();

} // _closeDataStep


// ------------------------------------------------------------------------------------------------
// Write output for step in solution.
void
pylith::meshio::OutputPhysics::_writeDataStep(const PylithReal t,
                                              const PylithInt tindex,
                                              const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputPhysics::_writeDataStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.getLabel()<<")");

    assert(_physics);
    const pylith::topology::Field* auxiliaryField = _physics->getAuxiliaryField();
    const pylith::topology::Field* derivedField = _physics->getDerivedField();
    const pylith::topology::Mesh& domainMesh = _physics->getPhysicsDomainMesh();

    const pylith::string_vector& dataNames = _expandDataFieldNames(solution, auxiliaryField, derivedField);

    _openDataStep(t, domainMesh);

    if (auxiliaryField) { auxiliaryField->scatterLocalToOutput(); }
    PetscVec auxiliaryVector = (auxiliaryField) ? auxiliaryField->getOutputVector() : NULL;

    if (derivedField) { derivedField->scatterLocalToOutput(); }
    PetscVec derivedVector = (derivedField) ? derivedField->getOutputVector() : NULL;

    PetscVec solutionVector = solution.getOutputVector();assert(solutionVector);

    const char* labelName = _physics->getPhysicsLabelName();
    const int labelValue = _physics->getPhysicsLabelValue();

    const size_t numDataFields = dataNames.size();
    for (size_t i = 0; i < numDataFields; i++) {
        OutputSubfield* subfield = NULL;
        if (solution.hasSubfield(dataNames[i].c_str())) {
            subfield = OutputObserver::_getSubfield(solution, domainMesh, dataNames[i].c_str());assert(subfield);
            subfield->setLabel(labelName, labelValue);
            subfield->projectWithLabel(solutionVector);
        } else if (auxiliaryField && auxiliaryField->hasSubfield(dataNames[i].c_str())) {
            subfield = OutputObserver::_getSubfield(*auxiliaryField, domainMesh, dataNames[i].c_str());assert(subfield);
            subfield->project(auxiliaryVector);
        } else if (derivedField && derivedField->hasSubfield(dataNames[i].c_str())) {
            subfield = OutputObserver::_getSubfield(*derivedField, domainMesh, dataNames[i].c_str());assert(subfield);
            subfield->setLabel(labelName, labelValue);
            subfield->project(derivedVector);
        } else {
            std::ostringstream msg;
            msg << "Internal Error: Could not find subfield '" << dataNames[i] << "' for data output.";
            PYLITH_COMPONENT_ERROR(msg.str());
            throw std::runtime_error(msg.str());
        } // if/else

        OutputObserver::_appendField(t, *subfield);
    } // for
    _closeDataStep();

    PYLITH_METHOD_END;
} // _writeDataStep


// ------------------------------------------------------------------------------------------------
// Names of information fields for output.
pylith::string_vector
pylith::meshio::OutputPhysics::_expandInfoFieldNames(const pylith::topology::Field* auxField) const {
    PYLITH_METHOD_BEGIN;

    if (auxField && (1 == _infoFieldNames.size()) && (std::string("all") == _infoFieldNames[0])) {
        PYLITH_METHOD_RETURN(auxField->getSubfieldNames());
    } // if

    PYLITH_METHOD_RETURN(_infoFieldNames);
} // _expandInfoFieldNames


// ------------------------------------------------------------------------------------------------
// Names of data fields for output.
pylith::string_vector
pylith::meshio::OutputPhysics::_expandDataFieldNames(const pylith::topology::Field& solution,
                                                     const pylith::topology::Field* auxField,
                                                     const pylith::topology::Field* derivedField) const {
    PYLITH_METHOD_BEGIN;

    if ((1 == _dataFieldNames.size()) && (std::string("all") == _dataFieldNames[0])) {
        pylith::string_vector dataNames = pylith::topology::FieldOps::getSubfieldNamesDomain(solution);

        if (auxField) {
            const pylith::string_vector& auxSubfields = auxField->getSubfieldNames();
            const size_t origSize = dataNames.size();
            dataNames.resize(origSize + auxSubfields.size());

            size_t numAdd = 0;
            for (size_t i = 0, iName = origSize; i < auxSubfields.size(); ++i) {
                const pylith::topology::Field::SubfieldInfo& info = auxField->getSubfieldInfo(auxSubfields[i].c_str());
                if (info.description.hasHistory) {
                    dataNames[iName++] = auxSubfields[i];
                    ++numAdd;
                } // if
            } // for
            dataNames.resize(origSize + numAdd);
        } // if
        if (derivedField) {
            const pylith::string_vector& derivedSubfields = derivedField->getSubfieldNames();
            const size_t origSize = dataNames.size();
            const size_t numAdd = derivedSubfields.size();
            dataNames.resize(origSize + numAdd);
            for (size_t iAdd = 0, iName = origSize; iAdd < numAdd; ++iAdd, ++iName) {
                dataNames[iName] = derivedSubfields[iAdd];
            } // for
        } // if
        PYLITH_METHOD_RETURN(dataNames);
    } // if

    PYLITH_METHOD_RETURN(_dataFieldNames);
} // _expandDataFieldNames


// End of file
