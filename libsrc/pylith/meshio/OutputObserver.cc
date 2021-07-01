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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "OutputObserver.hh" // Implementation of class methods

#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/meshio/OutputTrigger.hh" // USES OutputTrigger
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <iostream> // USES std::cout
#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputObserver::OutputObserver(void) :
    _timeScale(1.0),
    _writer(NULL),
    _trigger(NULL),
    _outputBasisOrder(1)
{}


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
// Get output subfield, creating if necessary.
pylith::meshio::OutputSubfield*
pylith::meshio::OutputObserver::_getSubfield(const pylith::topology::Field& field,
                                             const pylith::topology::Mesh& submesh,
                                             const char* name) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_getSubfield(field="<<field.getLabel()<<", name="<<name<<", submesh="<<typeid(submesh).name()<<")");

    if (_subfields.count(name) == 0) {
        _subfields[name] = OutputSubfield::create(field, submesh, name, _outputBasisOrder);
    } // if

    PYLITH_METHOD_RETURN(_subfields[name]);
} // _getSubfield


// ------------------------------------------------------------------------------------------------
// Append finite-element vertex field to file.
void
pylith::meshio::OutputObserver::_appendField(const PylithReal t,
                                             const pylith::meshio::OutputSubfield& subfield) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_appendField(t="<<t<<", subfield="<<typeid(subfield).name()<<")");

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

    PYLITH_METHOD_END;
} // _appendField


// End of file
