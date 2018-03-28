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

#include "OutputManager.hh" // Implementation of class methods

#include "DataWriter.hh" // USES DataWriter
#include "VertexFilter.hh" // USES VertexFilter
#include "CellFilter.hh" // USES CellFilter

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields

#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <iostream> // USES std::cout
#include <typeinfo> \
    // USES typeid()

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputManager::OutputManager(void) :
    _fields(NULL),
    _writer(NULL),
    _vertexFilter(NULL),
    _cellFilter(NULL),
    _timeSkip(0.0),
    _timeWrote(-PYLITH_MAXSCALAR),
    _numTimeStepsSkip(0),
    _timeStepWrote(PYLITH_MININT+10),
    _trigger(SKIP_TIMESTEPS)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputManager::~OutputManager(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputManager::deallocate(void) {
    _writer = NULL; // :TODO: Use shared pointer
    _vertexFilter = NULL; // :TODO: Use shared pointer
    _cellFilter = NULL; // :TODO: Use shared pointer
    delete _fields; _fields = NULL;
} // deallocate

// ----------------------------------------------------------------------
// Set trigger for how often to write output.
void
pylith::meshio::OutputManager::trigger(TriggerEnum flag) {
    PYLITH_COMPONENT_DEBUG("OutputManager::trigger(flag="<<flag<<")");

    _trigger = flag;
} // trigger

// ----------------------------------------------------------------------
// Get trigger for how often to write otuput.
pylith::meshio::OutputManager::TriggerEnum
pylith::meshio::OutputManager::trigger(void) const {
    return _trigger;
}  // trigger

// ----------------------------------------------------------------------
// Set number of time steps to skip between writes.
void
pylith::meshio::OutputManager::numTimeStepsSkip(const int value) {
    PYLITH_COMPONENT_DEBUG("OutputManager::numTimeStepsSkip(value="<<value<<")");

    _numTimeStepsSkip = (value >= 0) ? value : 0;
} // numTimeStepsSkip

// ----------------------------------------------------------------------
// Get number of time steps to skip between writes.
int
pylith::meshio::OutputManager::numTimeStepsSkip(void) const {
    return _numTimeStepsSkip;
} // numTimeStepsSkip

// ----------------------------------------------------------------------
// Set elapsed time between writes.
void
pylith::meshio::OutputManager::timeSkip(const double value) {
    PYLITH_COMPONENT_DEBUG("OutputManager::timeSkip(value="<<value<<")");

    _timeSkip = (value >= 0.0) ? value : 0.0;
} // timeSkip

// ----------------------------------------------------------------------
// Get elapsed time between writes.
double
pylith::meshio::OutputManager::timeSkip(void) const {
    return _timeSkip;
} // timeSkip

// ----------------------------------------------------------------------
// Set writer to write data to file.
void
pylith::meshio::OutputManager::writer(DataWriter* const datawriter) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::write(datawriter="<<typeid(datawriter).name()<<")");

    _writer = datawriter; // :TODO: Use shared pointer

    PYLITH_METHOD_END;
} // writer

// ----------------------------------------------------------------------
// Set filter for vertex data.
void
pylith::meshio::OutputManager::vertexFilter(VertexFilter* const filter) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::vertexFilter(filter="<<typeid(filter).name()<<")");

    _vertexFilter = filter; // :TODO: Use shared pointer

    PYLITH_METHOD_END;
} // vertexFilter

// ----------------------------------------------------------------------
// Set filter for cell data.
void
pylith::meshio::OutputManager::cellFilter(CellFilter* const filter) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::cellFilter(filter="<<typeid(filter).name()<<")");

    _cellFilter = filter; // :TODO: Use shared pointer

    PYLITH_METHOD_END;
} // cellFilter

// ----------------------------------------------------------------------
// Prepare for output.
void
pylith::meshio::OutputManager::open(const pylith::topology::Mesh& mesh,
                                    const bool isInfo,
                                    const char* label,
                                    const int labelId) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::open(mesh="<<typeid(mesh).name()<<", isInfo="<<isInfo<<", label="<<(label ? label : "NULL")<<", labelId="<<labelId<<")");


    if (!_writer) {
        std::ostringstream msg;
        if (label) {
            msg << "Writer for output manager for " << label << " not set.";
            throw std::runtime_error(msg.str());
        } else {
            throw std::runtime_error("Writer for output manager not set.");
        } // if/else
    } // if

    assert(_writer);
    _writer->open(mesh, isInfo, label, labelId);

    PYLITH_METHOD_END;
} // open

// ----------------------------------------------------------------------
/// Close output files.
void
pylith::meshio::OutputManager::close(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::close()");

    assert(_writer);
    _writer->close();

    PYLITH_METHOD_END;
} // close

// ----------------------------------------------------------------------
// Setup file for writing fields at time step.
void
pylith::meshio::OutputManager::openTimeStep(const PylithReal t,
                                            const pylith::topology::Mesh& mesh,
                                            const char* label,
                                            const int labelId) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::openTimeStep(t="<<t<<", mesh="<<typeid(mesh).name()<<", label="<<(label ? label : "NULL")<<", labelId="<<labelId<<")");

    assert(_writer);
    _writer->openTimeStep(t, mesh, label, labelId);

    PYLITH_METHOD_END;
} // openTimeStep

// ----------------------------------------------------------------------
// End writing fields at time step.
void
pylith::meshio::OutputManager::closeTimeStep(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::closeTimeStep()");

    assert(_writer);
    _writer->closeTimeStep();

    PYLITH_METHOD_END;
} // closeTimeStep

// ----------------------------------------------------------------------
// Append finite-element vertex field to file.
void
pylith::meshio::OutputManager::appendVertexField(const PylithReal t,
                                                 pylith::topology::Field& field,
                                                 const topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::appendVertexField(t="<<t<<", field="<<typeid(field).name()<<", mesh="<<typeid(mesh).name()<<")");

    topology::Field& fieldFiltered = (!_vertexFilter) ? field : _vertexFilter->filter(field);
    topology::Field& fieldDimensioned = _dimension(fieldFiltered);

    _writer->writeVertexField(t, fieldDimensioned, mesh);

    PYLITH_METHOD_END;
} // appendVertexField

// ----------------------------------------------------------------------
// Append finite-element cell field to file.
void
pylith::meshio::OutputManager::appendCellField(const PylithReal t,
                                               pylith::topology::Field& field,
                                               const char* label,
                                               const int labelId) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::appendCellField(t="<<t<<", field="<<typeid(field).name()<<", label="<<label<<", labelId="<<labelId<<")");

    topology::Field& fieldFiltered = (!_cellFilter) ? field : _cellFilter->filter(field, label, labelId);
    topology::Field& fieldDimensioned = _dimension(fieldFiltered);

    try {
        _writer->writeCellField(t, fieldDimensioned, label, labelId);
    } catch (std::runtime_error e) {
        std::cout << "ERROR: " << e.what() << std::endl<<std::endl<<std::endl;
    } // try/catch

    PYLITH_METHOD_END;
} // appendCellField

// ----------------------------------------------------------------------
// Check whether we want to write output at time t.
bool
pylith::meshio::OutputManager::shouldWrite(const PylithReal t,
                                           const PylithInt timeStep) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::shouldWrite(t="<<t<<", timeStep="<<timeStep<<")");

    bool shouldWrite = false;
    switch (_trigger) {
    case SKIP_TIMESTEPS:
        if (timeStep - _timeStepWrote > _numTimeStepsSkip) {
            shouldWrite = true;
            _timeStepWrote = timeStep;
        } // if
        break;
    case ELAPSED_TIME:
        if (t - _timeWrote >= _timeSkip) {
            shouldWrite = true;
            _timeWrote = t;
        } // if
        break;
    default:
        PYLITH_COMPONENT_ERROR("Unknown trigger type.");
        throw std::logic_error("Unknown trigger type.");
    } // switch

    PYLITH_METHOD_RETURN(shouldWrite);
} // shouldWrite

// ----------------------------------------------------------------------
/** Get buffer for field.
 *
 * Find the most appropriate buffer that matches field, reusing and reallocating as necessary.
 *
 * @param[in] fieldIn Input field.
 * @param[in] name Name of subfield (optional).
 * @returns Field to use as buffer for outputting field.
 */
pylith::topology::Field&
pylith::meshio::OutputManager::getBuffer(const pylith::topology::Field& fieldIn,
                                         const char* name) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::_getBuffer(fieldIn="<<fieldIn.label()<<")");

    pylith::topology::FieldBase::VectorFieldEnum fieldType = pylith::topology::FieldBase::MULTI_OTHER;
    if (name) {
        const pylith::topology::Field::SubfieldInfo& info = fieldIn.subfieldInfo(name);
        fieldType = info.description.vectorFieldType;
    } else {
        //fieldType = fieldIn.vectorFieldType();
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
        throw std::logic_error("Unknown field type in OutputManager::_getBuffer().");
    } // switch

    if (!_fields) {
        _fields = new topology::Fields(fieldIn.mesh()); assert(_fields);
    } // if

    if (!_fields->hasField(fieldName.c_str())) {
        _fields->add(fieldName.c_str(), fieldIn.label());
        topology::Field& fieldOut = _fields->get(fieldName.c_str());
        if (!name) {
            fieldOut.cloneSection(fieldIn);
        } // if/else
          //fieldOut.vectorFieldType(fieldIn.vectorFieldType());
          //fieldOut.scale(fieldIn.scale());
    } // if
    topology::Field& fieldOut = _fields->get(fieldName.c_str());
    if (name) {
        fieldOut.copySubfield(fieldIn, name);
    } else {
        fieldOut.copy(fieldIn);
    } // if/else
    fieldOut.dimensionalizeOkay(true);

    PYLITH_METHOD_RETURN(fieldOut);
} // getBuffer

// ----------------------------------------------------------------------
// Dimension field.
pylith::topology::Field&
pylith::meshio::OutputManager::_dimension(pylith::topology::Field& fieldIn) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::_dimension(fieldIn="<<fieldIn.label()<<")");

#if 0 // :TODO: FIX THIS.
    if (1.0 == fieldIn.scale()) {
        PYLITH_METHOD_RETURN(fieldIn);
    }
#endif

    if (fieldIn.dimensionalizeOkay()) {
        fieldIn.dimensionalize();
        PYLITH_METHOD_RETURN(fieldIn);
    } else {
        pylith::topology::Field& fieldOut = getBuffer(fieldIn);
        fieldOut.copy(fieldIn);
        fieldOut.dimensionalizeOkay(true);
        fieldOut.dimensionalize();

        PYLITH_METHOD_RETURN(fieldOut);
    } // if/else

    // Satisfy return value. Should never get this far.
    PYLITH_METHOD_RETURN(fieldIn);
} // _dimension


// End of file
