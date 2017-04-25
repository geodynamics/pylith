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

#include "OutputManagerNew.hh" // Implementation of class methods

#include "DataWriter.hh" // USES DataWriter
#include "VertexFilter.hh" // USES VertexFilter
#include "CellFilter.hh" // USES CellFilter

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields

#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <iostream> // USES std::cout
#include <typeinfo> // USES typeid()

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputManagerNew::OutputManagerNew(void) :
    _fields(NULL),
    _coordsys(NULL),
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
pylith::meshio::OutputManagerNew::~OutputManagerNew(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputManagerNew::deallocate(void)
{ // deallocate
    _writer = NULL; // :TODO: Use shared pointer
    _vertexFilter = NULL; // :TODO: Use shared pointer
    _cellFilter = NULL; // :TODO: Use shared pointer
    delete _coordsys; _coordsys = NULL;
    delete _fields; _fields = NULL;
} // deallocate

// ----------------------------------------------------------------------
// Set trigger for how often to write output.
void
pylith::meshio::OutputManagerNew::trigger(TriggerEnum flag)
{ // trigger
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::trigger(flag="<<flag<<")");

    _trigger = flag;
} // trigger

// ----------------------------------------------------------------------
// Get trigger for how often to write otuput.
pylith::meshio::OutputManagerNew::TriggerEnum
pylith::meshio::OutputManagerNew::trigger(void) const
{  // trigger
    return _trigger;
}  // trigger

// ----------------------------------------------------------------------
// Set number of time steps to skip between writes.
void
pylith::meshio::OutputManagerNew::numTimeStepsSkip(const int value)
{ // numTimeStepsSkip
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::numTimeStepsSkip(value="<<value<<")");

    _numTimeStepsSkip = (value >= 0) ? value : 0;
} // numTimeStepsSkip

// ----------------------------------------------------------------------
// Get number of time steps to skip between writes.
int
pylith::meshio::OutputManagerNew::numTimeStepsSkip(void) const
{ // numTimeStepsSkip
    return _numTimeStepsSkip;
} // numTimeStepsSkip

// ----------------------------------------------------------------------
// Set elapsed time between writes.
void
pylith::meshio::OutputManagerNew::timeSkip(const double value)
{ // timeSkip
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::timeSkip(value="<<value<<")");

    _timeSkip = (value >= 0.0) ? value : 0.0;
} // timeSkip

// ----------------------------------------------------------------------
// Get elapsed time between writes.
double
pylith::meshio::OutputManagerNew::timeSkip(void) const
{ // timeSkip
    return _timeSkip;
} // timeSkip

// ----------------------------------------------------------------------
// Set coordinate system in output. The vertex fields in the output
void
pylith::meshio::OutputManagerNew::coordsys(const spatialdata::geocoords::CoordSys* cs)
{ // coordsys
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::coordsys(cs="<<typeid(*cs).name()<<")");

    PYLITH_JOURNAL_WARNING("Coordinate system ignored in OutputManagerNew.");

    delete _coordsys; _coordsys = (cs) ? cs->clone() : 0;

    PYLITH_METHOD_END;
} // coordsys

// ----------------------------------------------------------------------
// Set writer to write data to file.
void
pylith::meshio::OutputManagerNew::writer(DataWriter* const datawriter)
{ // writer
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::write(datawriter="<<typeid(datawriter).name()<<")");

    _writer = datawriter; // :TODO: Use shared pointer

    PYLITH_METHOD_END;
} // writer

// ----------------------------------------------------------------------
// Set filter for vertex data.
void
pylith::meshio::OutputManagerNew::vertexFilter(VertexFilter* const filter)
{ // vertexFilter
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::vertexFilter(filter="<<typeid(filter).name()<<")");

    _vertexFilter = filter; // :TODO: Use shared pointer

    PYLITH_METHOD_END;
} // vertexFilter

// ----------------------------------------------------------------------
// Set filter for cell data.
void
pylith::meshio::OutputManagerNew::cellFilter(CellFilter* const filter)
{ // cellFilter
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::cellFilter(filter="<<typeid(filter).name()<<")");

    _cellFilter = filter; // :TODO: Use shared pointer

    PYLITH_METHOD_END;
} // cellFilter

// ----------------------------------------------------------------------
// Prepare for output.
void
pylith::meshio::OutputManagerNew::open(const topology::Mesh& mesh,
                                       const bool isInfo,
                                       const char* label,
                                       const int labelId)
{ // open
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::open(mesh="<<typeid(mesh).name()<<", isInfo="<<isInfo<<", label="<<(label ? label : "NULL")<<", labelId="<<labelId<<")");


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
pylith::meshio::OutputManagerNew::close(void)
{ // close
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::close()");

    assert(_writer);
    _writer->close();

    PYLITH_METHOD_END;
} // close

// ----------------------------------------------------------------------
// Setup file for writing fields at time step.
void
pylith::meshio::OutputManagerNew::openTimeStep(const PylithReal t,
                                               const topology::Mesh& mesh,
                                               const char* label,
                                               const int labelId)
{ // openTimeStep
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::openTimeStep(t="<<t<<", mesh="<<typeid(mesh).name()<<", label="<<(label ? label : "NULL")<<", labelId="<<labelId<<")");

    assert(_writer);
    _writer->openTimeStep(t, mesh, label, labelId);

    PYLITH_METHOD_END;
} // openTimeStep

// ----------------------------------------------------------------------
// End writing fields at time step.
void
pylith::meshio::OutputManagerNew::closeTimeStep(void)
{ // closeTimeStep
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::closeTimeStep()");

    assert(_writer);
    _writer->closeTimeStep();

    PYLITH_METHOD_END;
} // closeTimeStep

// ----------------------------------------------------------------------
// Append finite-element vertex field to file.
void
pylith::meshio::OutputManagerNew::appendVertexField(const PylithReal t,
                                                    topology::Field& field,
                                                    const topology::Mesh& mesh)
{ // appendVertexField
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::appendVertexField(t="<<t<<", field="<<typeid(field).name()<<", mesh="<<typeid(mesh).name()<<")");

    topology::Field& fieldFiltered = (!_vertexFilter) ? field : _vertexFilter->filter(field);
    topology::Field& fieldDimensioned = _dimension(fieldFiltered);

    _writer->writeVertexField(t, fieldDimensioned, mesh);

    PYLITH_METHOD_END;
} // appendVertexField

// ----------------------------------------------------------------------
// Append finite-element cell field to file.
void
pylith::meshio::OutputManagerNew::appendCellField(const PylithReal t,
                                                  topology::Field& field,
                                                  const char* label,
                                                  const int labelId)
{ // appendCellField
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::appendCellField(t="<<t<<", field="<<typeid(field).name()<<", label="<<label<<", labelId="<<labelId<<")");

    topology::Field& fieldFiltered = (!_cellFilter) ? field : _cellFilter->filter(field, label, labelId);
    topology::Field& fieldDimensioned = _dimension(fieldFiltered);

    try {
        _writer->writeCellField(t, fieldDimensioned, label, labelId);
    } catch(std::runtime_error e) {
        std::cout << "ERROR: " << e.what() << std::endl<<std::endl<<std::endl;
    } // try/catch

    PYLITH_METHOD_END;
} // appendCellField

// ----------------------------------------------------------------------
// Check whether we want to write output at time t.
bool
pylith::meshio::OutputManagerNew::shouldWrite(const PylithReal t,
                                              const PylithInt timeStep)
{ // shouldWrite
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::shouldWrite(t="<<t<<", timeStep="<<timeStep<<")");

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
        PYLITH_JOURNAL_ERROR("Unknown trigger type.");
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
pylith::meshio::OutputManagerNew::getBuffer(const pylith::topology::Field& fieldIn,
                                            const char* name)
{ // getBuffer
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::_getBuffer(fieldIn="<<fieldIn.label()<<")");

    pylith::topology::FieldBase::VectorFieldEnum fieldType = pylith::topology::FieldBase::MULTI_OTHER;
    if (name) {
        const pylith::topology::Field::SubfieldInfo& info = fieldIn.subfieldInfo(name);
        fieldType = info.metadata.vectorFieldType;
    } else {
        fieldType = fieldIn.vectorFieldType();
    } // if/else

    std::string fieldName = "buffer (other)";
    switch (fieldType)
    { // switch
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
        PYLITH_JOURNAL_ERROR("Unknown field type '"<<fieldIn.vectorFieldType()<<"' for field '"<<fieldIn.label()<<"'.");
        throw std::logic_error("Unknown field type in OutputManagerNew::_getBuffer().");
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
        fieldOut.vectorFieldType(fieldIn.vectorFieldType());
        fieldOut.scale(fieldIn.scale());
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
pylith::meshio::OutputManagerNew::_dimension(topology::Field& fieldIn)
{ // _dimension
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("OutputManagerNew::_dimension(fieldIn="<<fieldIn.label()<<")");

    if (1.0 == fieldIn.scale())
        PYLITH_METHOD_RETURN(fieldIn);

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
