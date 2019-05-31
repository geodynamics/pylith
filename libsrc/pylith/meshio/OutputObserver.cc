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

#include "OutputObserver.hh" // Implementation of class methods

#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/meshio/FieldFilter.hh" // USES FieldFilter
#include "pylith/meshio/OutputTrigger.hh" // USES OutputTrigger

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields

#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <iostream> // USES std::cout
#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputObserver::OutputObserver(void) :
    _timeScale(1.0),
    _fields(NULL),
    _writer(NULL),
    _fieldFilter(NULL),
    _trigger(NULL)
{}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputObserver::~OutputObserver(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputObserver::deallocate(void) {
    if (_writer) {
        _writer->close();
        _writer->deallocate();
    }
    if (_fieldFilter) { _fieldFilter->deallocate(); }

    _writer = NULL; // :TODO: Use shared pointer
    _fieldFilter = NULL; // :TODO: Use shared pointer
    _trigger = NULL; // :TODO: Use shared pointer

    delete _fields;_fields = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set trigger for how often to write output.
void
pylith::meshio::OutputObserver::setTrigger(pylith::meshio::OutputTrigger* const trigger) {
    PYLITH_COMPONENT_DEBUG("OutputObserver::setTrigger(otrigger="<<typeid(trigger).name()<<")");

    _trigger = trigger;
} // setTrigger


// ---------------------------------------------------------------------------------------------------------------------
// Get trigger for how often to write otuput.
const pylith::meshio::OutputTrigger*
pylith::meshio::OutputObserver::getTrigger(void) const {
    return _trigger;
} // getTrigger


// ---------------------------------------------------------------------------------------------------------------------
// Set writer to write data to file.
void
pylith::meshio::OutputObserver::setWriter(DataWriter* const writer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputObserver::setWrite(datawriter="<<typeid(writer).name()<<")");

    _writer = writer; // :TODO: Use shared pointer

    PYLITH_METHOD_END;
} // setWriter


// ---------------------------------------------------------------------------------------------------------------------
// Set filter for fields.
void
pylith::meshio::OutputObserver::setFieldFilter(FieldFilter* const filter) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputObserver::setFieldFilter(filter="<<typeid(filter).name()<<")");

    _fieldFilter = filter; // :TODO: Use shared pointer

    PYLITH_METHOD_END;
} // setFieldFilter


// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
/** Get buffer for field.
 *
 * Find the most appropriate buffer that matches field, reusing and reallocating as necessary.
 *
 * @param[in] fieldIn Input field.
 * @param[in] name Name of subfield (optional).
 * @returns Field to use as buffer for outputting field.
 */
pylith::topology::Field*
pylith::meshio::OutputObserver::_getBuffer(const pylith::topology::Field& fieldIn,
                                           const char* name) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputObserver::_getBuffer(fieldIn="<<fieldIn.label()<<")");

    pylith::topology::FieldBase::VectorFieldEnum fieldType = pylith::topology::FieldBase::MULTI_OTHER;
    if (name) {
        fieldType = fieldIn.subfieldInfo(name).description.vectorFieldType;
    } else {
        // Get vector field type for subfield if only one subfield in field.
        const pylith::string_vector& subfieldNames = fieldIn.subfieldNames();
        if (size_t(1) == subfieldNames.size()) {
            fieldType = fieldIn.subfieldInfo(subfieldNames[0].c_str()).description.vectorFieldType;
        } else {
            PYLITH_COMPONENT_ERROR("No subfield specified for field '"<<fieldIn.label() <<"' with multiple subfields.");
            throw std::runtime_error("No subfield specified for field with multiple fields.");
        } // if/else
    } // if/else

    std::string fieldName = "buffer (other)";
    switch (fieldType) { // switch
    case topology::FieldBase::SCALAR:
        fieldName = "buffer (scalar)";
        break;
    case topology::FieldBase::VECTOR:
        fieldName = "buffer (vector)";
        break;
    case topology::FieldBase::TENSOR:
        fieldName = "buffer (tensor)";
        break;
    case topology::FieldBase::OTHER:
        fieldName = "buffer (other)";
        break;
    case topology::FieldBase::MULTI_SCALAR:
        fieldName = "buffer (multiple scalars)";
        break;
    case topology::FieldBase::MULTI_VECTOR:
        fieldName = "buffer (multiple vectors)";
        break;
    case topology::FieldBase::MULTI_TENSOR:
        fieldName = "buffer (multiple tensors)";
        break;
    case topology::FieldBase::MULTI_OTHER:
        fieldName = "buffer (multiple others)";
        break;
    default:
        PYLITH_COMPONENT_ERROR("Unknown field type '"<<fieldType<<"' for field '"<<fieldIn.label()<<"'.");
        throw std::logic_error("Unknown field type in OutputObserver::_getBuffer().");
    } // switch

    delete _fields;_fields = NULL; // :KLUDGE: :TODO: @brad DS is not getting set when extracting subfield.
    if (!_fields) {
        _fields = new topology::Fields(fieldIn.mesh());assert(_fields);
    } // if

    if (!_fields->hasField(fieldName.c_str())) {
        _fields->add(fieldName.c_str(), fieldIn.label());
        topology::Field& fieldOut = _fields->get(fieldName.c_str());
        if (!name) {
            fieldOut.cloneSection(fieldIn);
        } // if/else
          // fieldOut.vectorFieldType(fieldIn.vectorFieldType());
          // fieldOut.scale(fieldIn.scale());
    } // if
    pylith::topology::Field& fieldOut = _fields->get(fieldName.c_str());
    if (name) {
        fieldOut.copySubfield(fieldIn, name);
    } else {
        fieldOut.copy(fieldIn);
    } // if/else
    fieldOut.dimensionalizeOkay(true);

    PYLITH_METHOD_RETURN(&fieldOut);
} // _getBuffer


// ---------------------------------------------------------------------------------------------------------------------
// Dimension field.
pylith::topology::Field*
pylith::meshio::OutputObserver::_dimensionField(pylith::topology::Field* fieldIn) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputObserver::_dimension(fieldIn="<<typeid(fieldIn).name()<<")");

    if (!fieldIn) { PYLITH_METHOD_RETURN(NULL); }

    assert(fieldIn);

    // Check to see if all subfields have scales of 1.0.
    bool needDimensioning = false;
    const pylith::string_vector& subfieldNames = fieldIn->subfieldNames();
    const size_t numSubfields = subfieldNames.size();
    for (size_t i = 0; i < numSubfields; ++i) {
        if (fieldIn->subfieldInfo(subfieldNames[i].c_str()).description.scale != 1.0) {
            needDimensioning = true;
            break;
        } // if
    } // for
    if (!needDimensioning) { PYLITH_METHOD_RETURN(fieldIn); }

    if (fieldIn->dimensionalizeOkay()) {
        fieldIn->dimensionalize();
        PYLITH_METHOD_RETURN(fieldIn);
    } else {
        pylith::topology::Field* fieldOut = _getBuffer(*fieldIn);
        fieldOut->copy(*fieldIn);
        fieldOut->dimensionalizeOkay(true);
        fieldOut->dimensionalize();

        PYLITH_METHOD_RETURN(fieldOut);
    } // if/else

    // Satisfy return value. Should never get this far.
    PYLITH_METHOD_RETURN(fieldIn);
} // _dimensionField


// ---------------------------------------------------------------------------------------------------------------------
// Get basis order of field.
int
pylith::meshio::OutputObserver::_getBasisOrder(const pylith::topology::Field& field) {
    PYLITH_METHOD_BEGIN;

    int basisOrder = -1;

    const pylith::string_vector& subfieldNames = field.subfieldNames();
    const size_t numSubfields = subfieldNames.size();
    if (1 == numSubfields) {
        basisOrder = field.subfieldInfo(subfieldNames[0].c_str()).fe.basisOrder;
    } else {
        PYLITH_COMPONENT_ERROR("Expected one subfield in field '"<<field.label()<<"'.");
    } // if/else

    PYLITH_METHOD_RETURN(basisOrder);
} // _getBasisOrder


// End of file
