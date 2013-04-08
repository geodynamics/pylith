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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "DataWriter.hh" // USES DataWriter
#include "VertexFilter.hh" // USES VertexFilter
#include "CellFilter.hh" // USES CellFilter

#include "pylith/topology/Fields.hh" // USES Fields

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type, typename field_type>
pylith::meshio::OutputManager<mesh_type, field_type>::OutputManager(void) :
  _coordsys(0),
  _writer(0),
  _vertexFilter(0),
  _cellFilter(0),
  _fields(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type, typename field_type>
pylith::meshio::OutputManager<mesh_type, field_type>::~OutputManager(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::deallocate(void)
{ // deallocate
  _writer = 0; // :TODO: Use shared pointer
  _vertexFilter = 0; // :TODO: Use shared pointer
  _cellFilter = 0; // :TODO: Use shared pointer
  delete _coordsys; _coordsys = 0;
  delete _fields; _fields = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Set coordinate system in output. The vertex fields in the output
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::coordsys(const spatialdata::geocoords::CoordSys* cs)
{ // coordsys
  PYLITH_METHOD_BEGIN;

  delete _coordsys; _coordsys = (cs) ? cs->clone() : 0;

  PYLITH_METHOD_END;
} // coordsys

// ----------------------------------------------------------------------
// Set writer to write data to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::writer(DataWriter<mesh_type, field_type>* const datawriter)
{ // writer
  PYLITH_METHOD_BEGIN;

  _writer = datawriter; // :TODO: Use shared pointer

  PYLITH_METHOD_END;
} // writer

// ----------------------------------------------------------------------
// Set filter for vertex data.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::vertexFilter(VertexFilter<field_type>* const filter)
{ // vertexFilter
  PYLITH_METHOD_BEGIN;

  _vertexFilter = filter; // :TODO: Use shared pointer

  PYLITH_METHOD_END;
} // vertexFilter

// ----------------------------------------------------------------------
// Set filter for cell data.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::cellFilter(CellFilter<mesh_type, field_type>* const filter)
{ // cellFilter
  PYLITH_METHOD_BEGIN;

  _cellFilter = filter; // :TODO: Use shared pointer

  PYLITH_METHOD_END;
} // cellFilter

// ----------------------------------------------------------------------
// Get fields used in output.
template<typename mesh_type, typename field_type>
const pylith::topology::Fields<field_type>*
pylith::meshio::OutputManager<mesh_type, field_type>::fields(void) const
{ // fields
  PYLITH_METHOD_BEGIN;

  PYLITH_METHOD_RETURN(_fields);
} // fields

// ----------------------------------------------------------------------
// Prepare for output.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::open(const mesh_type& mesh,
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
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::close(void)
{ // close
  PYLITH_METHOD_BEGIN;

  assert(_writer);
  _writer->close();

  PYLITH_METHOD_END;
} // close

// ----------------------------------------------------------------------
// Setup file for writing fields at time step.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::openTimeStep(const PylithScalar t,
								   const mesh_type& mesh,
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
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::closeTimeStep(void)
{ // closeTimeStep
  PYLITH_METHOD_BEGIN;

  assert(_writer);
  _writer->closeTimeStep();

  PYLITH_METHOD_END;
} // closeTimeStep

// ----------------------------------------------------------------------
// Append finite-element vertex field to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::appendVertexField(const PylithScalar t,
									field_type& field,
									const mesh_type& mesh)
{ // appendVertexField
  PYLITH_METHOD_BEGIN;

  field_type& fieldFiltered = (!_vertexFilter) ? field : _vertexFilter->filter(field);
  field_type& fieldDimensioned = _dimension(fieldFiltered);
  
  _writer->writeVertexField(t, fieldDimensioned, mesh);

  PYLITH_METHOD_END;
} // appendVertexField

// ----------------------------------------------------------------------
// Append finite-element cell field to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::appendCellField(const PylithScalar t,
								      field_type& field,
								      const char* label,
								      const int labelId)
{ // appendCellField
  PYLITH_METHOD_BEGIN;

  field_type& fieldFiltered = (!_cellFilter) ? field : _cellFilter->filter(field, label, labelId);
  field_type& fieldDimensioned = _dimension(fieldFiltered);

  try {
    _writer->writeCellField(t, fieldDimensioned, label, labelId);
  } catch(std::runtime_error e) {
    std::cout << "ERROR: " << e.what() << std::endl<<std::endl<<std::endl;
  } // try/catch

  PYLITH_METHOD_END;
} // appendCellField

// ----------------------------------------------------------------------
// Dimension field.
template<typename mesh_type, typename field_type>
field_type&
pylith::meshio::OutputManager<mesh_type, field_type>::_dimension(field_type& fieldIn)
{ // _dimension
  PYLITH_METHOD_BEGIN;

  if (1.0 == fieldIn.scale())
    PYLITH_METHOD_RETURN(fieldIn);

  if (fieldIn.addDimensionOkay()) {
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
	// Spit out usefule error message and stop via assert. If
	// optimized, throw exception.
	std::cerr << "Unknown field type '" << fieldIn.vectorFieldType()
		  << "'";
	assert(0);
	throw std::logic_error("Unknown field type");
      } // switch
    
    if (!_fields)
      _fields = new topology::Fields<field_type>(fieldIn.mesh());
    
    if (!_fields->hasField(fieldName.c_str())) {
      _fields->add(fieldName.c_str(), fieldIn.label());
      field_type& fieldOut = _fields->get(fieldName.c_str());
      fieldOut.cloneSection(fieldIn);
      fieldOut.vectorFieldType(fieldIn.vectorFieldType());
      fieldOut.scale(fieldIn.scale());
    } // if
    field_type& fieldOut = _fields->get(fieldName.c_str());
    fieldOut.copy(fieldIn);
    fieldOut.addDimensionOkay(true);
    fieldOut.dimensionalize();

    PYLITH_METHOD_RETURN(fieldOut);
  } // if/else

  // Satisfy return value. Should never get this far.
  PYLITH_METHOD_RETURN(fieldIn);
} // _dimension


// End of file 
