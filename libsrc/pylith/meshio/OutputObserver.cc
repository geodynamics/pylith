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

#include "pylith/meshio/OutputObserver.hh" // Implementation of class methods

#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/meshio/OutputTrigger.hh" // USES OutputTrigger
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/utils/constants.hh" // USES pylith::max_real
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include <iostream> // USES std::cout
#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace meshio {
        class _OutputObserver {
public:

            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static PylithInt getSubfield;
                static PylithInt appendField;
            };

        }; // _OutputObserver
    } // meshio
} // pylith

pylith::utils::EventLogger pylith::meshio::_OutputObserver::Events::logger;
PylithInt pylith::meshio::_OutputObserver::Events::getSubfield;
PylithInt pylith::meshio::_OutputObserver::Events::appendField;

// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_OutputObserver::Events::init(void) {
    logger.setClassName("OutputObserver");
    logger.initialize();
    getSubfield = logger.registerEvent("PL:OutputObserver:getSubfield");
    appendField = logger.registerEvent("PL:OutputObserver:appendField");
}


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputObserver::OutputObserver(void) :
    _timeScale(1.0),
    _outputMesh(NULL),
    _writer(NULL),
    _trigger(NULL),
    _outputBasisOrder(1) {
    _OutputObserver::Events::init();
}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputObserver::~OutputObserver(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputObserver::deallocate(void) {
    if (_writer) {
        _writer->close();
        _writer->deallocate();
    }

    typedef std::map<std::string, OutputSubfield*> subfield_t;
    for (subfield_t::iterator iter = _subfields.begin(); iter != _subfields.end(); ++iter) {
        delete iter->second;iter->second = NULL;
    } // for
    _subfields.clear();
    delete _outputMesh;_outputMesh = NULL;

    _writer = NULL; // :TODO: Use shared pointer
    _trigger = NULL; // :TODO: Use shared pointer

} // deallocate


// ------------------------------------------------------------------------------------------------
// Set trigger for how often to write output.
void
pylith::meshio::OutputObserver::setTrigger(pylith::meshio::OutputTrigger* const trigger) {
    PYLITH_COMPONENT_DEBUG("OutputObserver::setTrigger(otrigger="<<typeid(trigger).name()<<")");

    _trigger = trigger;
} // setTrigger


// ------------------------------------------------------------------------------------------------
// Get trigger for how often to write otuput.
const pylith::meshio::OutputTrigger*
pylith::meshio::OutputObserver::getTrigger(void) const {
    return _trigger;
} // getTrigger


// ------------------------------------------------------------------------------------------------
// Set writer to write data to file.
void
pylith::meshio::OutputObserver::setWriter(DataWriter* const writer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputObserver::setWrite(datawriter="<<typeid(writer).name()<<")");

    _writer = writer; // :TODO: Use shared pointer

    PYLITH_METHOD_END;
} // setWriter


// ------------------------------------------------------------------------------------------------
// Set basis order for output.
void
pylith::meshio::OutputObserver::setOutputBasisOrder(const int value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputObserver::setBasisOrder(value="<<value<<")");

    if (value < 0) {
        std::ostringstream msg;
        msg << "Basis order for output (" << value << ") must be nonnegative.";
        throw std::out_of_range(msg.str());
    } // if

    _outputBasisOrder = value;

    PYLITH_METHOD_END;
} // setOutputBasisOrder


// ------------------------------------------------------------------------------------------------

// Set number of refinement levels for output.
void
pylith::meshio::OutputObserver::setRefineLevels(const int value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputObserver::setRefineLevels(value="<<value<<")");

    if (value < 0) {
        std::ostringstream msg;
        msg << "Number of refinement levels for output (" << value << ") must be nonnegative.";
        throw std::out_of_range(msg.str());
    } // if

    _refineLevels = value;

    PYLITH_METHOD_END;
} // setRefineLevels


// ------------------------------------------------------------------------------------------------
// Set time scale.
void
pylith::meshio::OutputObserver::setTimeScale(const PylithReal value) {
    if (value <= 0.0) {
        std::ostringstream msg;
        msg << "Time scale ("<<value<<") for output observer is nonpositive.";
        throw std::logic_error(msg.str());
    } // if
    _timeScale = value;
} // setTimeScale


// ------------------------------------------------------------------------------------------------
// Get mesh associated with subfield output.
pylith::topology::Mesh*
pylith::meshio::OutputObserver::_getOutputMesh(const pylith::meshio::OutputSubfield& subfield) {
    PYLITH_METHOD_BEGIN;

    if (!_outputMesh) {
        _outputMesh = new pylith::topology::Mesh();
        PetscDM dmOutput = subfield.getOutputDM();
        PetscObjectReference((PetscObject) dmOutput);
        _outputMesh->setDM(dmOutput);
    } // if

    PYLITH_METHOD_RETURN(_outputMesh);
} // _getOutputMesh


// ------------------------------------------------------------------------------------------------
// Get output subfield, creating if necessary.
pylith::meshio::OutputSubfield*
pylith::meshio::OutputObserver::_getSubfield(const pylith::topology::Field& field,
                                             const pylith::topology::Mesh& submesh,
                                             const char* name) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_getSubfield(field="<<field.getLabel()<<", name="<<name<<", submesh="<<typeid(submesh).name()<<")");
    _OutputObserver::Events::logger.eventBegin(_OutputObserver::Events::getSubfield);

    if (0 == _subfields.count(name) ) {
        _subfields[name] = OutputSubfield::create(field, submesh, name, _outputBasisOrder, _refineLevels);
    } // if

    _OutputObserver::Events::logger.eventEnd(_OutputObserver::Events::getSubfield);
    PYLITH_METHOD_RETURN(_subfields[name]);
} // _getSubfield


// ------------------------------------------------------------------------------------------------
// Append finite-element vertex field to file.
void
pylith::meshio::OutputObserver::_appendField(const PylithReal t,
                                             const pylith::meshio::OutputSubfield& subfield) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_appendField(t="<<t<<", subfield="<<typeid(subfield).name()<<")");
    _OutputObserver::Events::logger.eventBegin(_OutputObserver::Events::appendField);

    // Use basis order from subfield since requested basis order for output may be greater than original basis order.
    const int basisOrder = subfield.getBasisOrder();
    switch (basisOrder) {
    case 0:
        _writer->writeCellField(t, subfield);
        break;

    case 1:
        _writer->writeVertexField(t, subfield);
        break;

    default:
        PYLITH_COMPONENT_ERROR(
            "Unsupported basis order ("<< basisOrder <<") for output. Skipping output of '"
                                       << subfield.getDescription().label << "' field."
            );
    } // switch

    _OutputObserver::Events::logger.eventEnd(_OutputObserver::Events::appendField);
    PYLITH_METHOD_END;
} // _appendField


// End of file
