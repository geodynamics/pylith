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
// Copyright (c) 2010-2015 University of California, Davis
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

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include <iostream> // USES std::cout

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputManager::OutputManager(void) :
  _fields(0),
  _coordsys(0),
  _writer(0),
  _vertexFilter(0),
  _cellFilter(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputManager::~OutputManager(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputManager::deallocate(void)
{ // deallocate
  _writer = 0; // :TODO: Use shared pointer
  _vertexFilter = 0; // :TODO: Use shared pointer
  _cellFilter = 0; // :TODO: Use shared pointer
  delete _coordsys; _coordsys = 0;
  delete _fields; _fields = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Set coordinate system in output. The vertex fields in the output
void
pylith::meshio::OutputManager::coordsys(const spatialdata::geocoords::CoordSys* cs)
{ // coordsys
  PYLITH_METHOD_BEGIN;

  delete _coordsys; _coordsys = (cs) ? cs->clone() : 0;

  PYLITH_METHOD_END;
} // coordsys

// ----------------------------------------------------------------------
// Set writer to write data to file.
void
pylith::meshio::OutputManager::writer(DataWriter* const datawriter)
{ // writer
  PYLITH_METHOD_BEGIN;

  _writer = datawriter; // :TODO: Use shared pointer

  PYLITH_METHOD_END;
} // writer

// ----------------------------------------------------------------------
// Set filter for vertex data.
void
pylith::meshio::OutputManager::vertexFilter(VertexFilter* const filter)
{ // vertexFilter
  PYLITH_METHOD_BEGIN;

  _vertexFilter = filter; // :TODO: Use shared pointer

  PYLITH_METHOD_END;
} // vertexFilter

// ----------------------------------------------------------------------
// Set filter for cell data.
void
pylith::meshio::OutputManager::cellFilter(CellFilter* const filter)
{ // cellFilter
  PYLITH_METHOD_BEGIN;

  _cellFilter = filter; // :TODO: Use shared pointer

  PYLITH_METHOD_END;
} // cellFilter

// ----------------------------------------------------------------------
// Get fields used in output.
const pylith::topology::Fields*
pylith::meshio::OutputManager::fields(void) const
{ // fields
  PYLITH_METHOD_BEGIN;

  PYLITH_METHOD_RETURN(_fields);
} // fields

// ----------------------------------------------------------------------
// Prepare for output.
void
pylith::meshio::OutputManager::open(const topology::Mesh& mesh,
				    const int numTimeSteps,
				    const char* label,
				    const int labelId)
{ // open
  PYLITH_METHOD_BEGIN;

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
  _writer->open(mesh, numTimeSteps, label, labelId);

  PYLITH_METHOD_END;
} // open

// ----------------------------------------------------------------------
/// Close output files.
void
pylith::meshio::OutputManager::close(void)
{ // close
  PYLITH_METHOD_BEGIN;

  assert(_writer);
  _writer->close();

  PYLITH_METHOD_END;
} // close

// ----------------------------------------------------------------------
// Setup file for writing fields at time step.
void
pylith::meshio::OutputManager::openTimeStep(const PylithScalar t,
					    const topology::Mesh& mesh,
					    const char* label,
					    const int labelId)
{ // openTimeStep
  PYLITH_METHOD_BEGIN;

  assert(_writer);
  _writer->openTimeStep(t, mesh, label, labelId);

  PYLITH_METHOD_END;
} // openTimeStep

// ----------------------------------------------------------------------
// End writing fields at time step.
void
pylith::meshio::OutputManager::closeTimeStep(void)
{ // closeTimeStep
  PYLITH_METHOD_BEGIN;

  assert(_writer);
  _writer->closeTimeStep();

  PYLITH_METHOD_END;
} // closeTimeStep

// ----------------------------------------------------------------------
// Append finite-element vertex field to file.
void
pylith::meshio::OutputManager::appendVertexField(const PylithScalar t,
						 topology::Field& field,
						 const topology::Mesh& mesh)
{ // appendVertexField
  PYLITH_METHOD_BEGIN;

  topology::Field& fieldFiltered = (!_vertexFilter) ? field : _vertexFilter->filter(field);
  topology::Field& fieldDimensioned = _dimension(fieldFiltered);
  
  _writer->writeVertexField(t, fieldDimensioned, mesh);

  PYLITH_METHOD_END;
} // appendVertexField

// ----------------------------------------------------------------------
// Append finite-element cell field to file.
void
pylith::meshio::OutputManager::appendCellField(const PylithScalar t,
					       topology::Field& field,
					       const char* label,
					       const int labelId)
{ // appendCellField
  PYLITH_METHOD_BEGIN;

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
// Dimension field.
pylith::topology::Field&
pylith::meshio::OutputManager::_dimension(topology::Field& fieldIn)
{ // _dimension
  PYLITH_METHOD_BEGIN;

  if (1.0 == fieldIn.scale())
    PYLITH_METHOD_RETURN(fieldIn);

  if (fieldIn.dimensionalizeOkay()) {
    fieldIn.dimensionalize();
    PYLITH_METHOD_RETURN(fieldIn);
  } else {
    std::string fieldName = "buffer (other)";
    switch (fieldIn.vectorFieldType())
      { // switch
      case topology::FieldBase::SCALAR :
	fieldName = "buffer (scalar)";
	break;
      case topology::FieldBase::VECTOR :
	fieldName = "buffer (vector)";
	break;
      case topology::FieldBase::TENSOR :
	fieldName = "buffer (tensor)";
	break;
      case topology::FieldBase::OTHER :
	fieldName = "buffer (other)";
	break;
      case topology::FieldBase::MULTI_SCALAR :
	fieldName = "buffer (multiple scalars)";
	break;
      case topology::FieldBase::MULTI_VECTOR :
	fieldName = "buffer (multiple vectors)";
	break;
      case topology::FieldBase::MULTI_TENSOR :
	fieldName = "buffer (multiple tensors)";
	break;
      case topology::FieldBase::MULTI_OTHER :
	fieldName = "buffer (multiple others)";
	break;
      default :
        // Spit out useful error message and stop via assert. If
        // optimized, throw exception.
        std::ostringstream msg;
        msg << "Unknown field type '" << fieldIn.vectorFieldType() << "'";
        throw std::logic_error(msg.str());
      } // switch
    
    if (!_fields) {
      _fields = new topology::Fields(fieldIn.mesh());assert(_fields);
    } // if
    
    if (!_fields->hasField(fieldName.c_str())) {
      _fields->add(fieldName.c_str(), fieldIn.label());
      topology::Field& fieldOut = _fields->get(fieldName.c_str());
      fieldOut.cloneSection(fieldIn);
      fieldOut.vectorFieldType(fieldIn.vectorFieldType());
      fieldOut.scale(fieldIn.scale());
    } // if
    topology::Field& fieldOut = _fields->get(fieldName.c_str());
    fieldOut.copy(fieldIn);
    fieldOut.dimensionalizeOkay(true);
    fieldOut.dimensionalize();

    PYLITH_METHOD_RETURN(fieldOut);
  } // if/else

  // Satisfy return value. Should never get this far.
  PYLITH_METHOD_RETURN(fieldIn);
} // _dimension


// End of file 
