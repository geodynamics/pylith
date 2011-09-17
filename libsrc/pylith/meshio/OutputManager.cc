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
// Copyright (c) 2010-2011 University of California, Davis
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
pylith::meshio::OutputManager<mesh_type, field_type>::coordsys(
				  const spatialdata::geocoords::CoordSys* cs)
{ // coordsys
  delete _coordsys; _coordsys = (0 != cs) ? cs->clone() : 0;
} // coordsys

// ----------------------------------------------------------------------
// Set writer to write data to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::writer(
			 DataWriter<mesh_type, field_type>* const datawriter)
{ // writer
  _writer = datawriter; // :TODO: Use shared pointer
} // writer

// ----------------------------------------------------------------------
// Set filter for vertex data.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::vertexFilter(
					VertexFilter<field_type>* const filter)
{ // vertexFilter
  _vertexFilter = filter; // :TODO: Use shared pointer
} // vertexFilter

// ----------------------------------------------------------------------
// Set filter for cell data.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::cellFilter(
			     CellFilter<mesh_type, field_type>* const filter)
{ // cellFilter
  _cellFilter = filter; // :TODO: Use shared pointer
} // cellFilter

// ----------------------------------------------------------------------
// Get fields used in output.
template<typename mesh_type, typename field_type>
const pylith::topology::Fields<field_type>*
pylith::meshio::OutputManager<mesh_type, field_type>::fields(void) const
{ // fields
  return _fields;
} // fields

// ----------------------------------------------------------------------
// Prepare for output.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::open(
						   const mesh_type& mesh,
						   const int numTimeSteps,
						   const char* label,
						   const int labelId)
{ // open
  assert(0 != _writer);

  _writer->open(mesh, numTimeSteps, label, labelId);
} // open

// ----------------------------------------------------------------------
/// Close output files.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::close(void)
{ // close
  assert(0 != _writer);
  _writer->close();
} // close

// ----------------------------------------------------------------------
// Setup file for writing fields at time step.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::openTimeStep(
						       const PylithScalar t,
						       const mesh_type& mesh,
						       const char* label,
						       const int labelId)
{ // openTimeStep
  assert(0 != _writer);
  _writer->openTimeStep(t, mesh, label, labelId);
} // openTimeStep

// ----------------------------------------------------------------------
// End writing fields at time step.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::closeTimeStep(void)
{ // closeTimeStep
  assert(0 != _writer);
  _writer->closeTimeStep();
} // closeTimeStep

// ----------------------------------------------------------------------
// Append finite-element vertex field to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::appendVertexField(
			                                const PylithScalar t,
							field_type& field,
							const mesh_type& mesh)
{ // appendVertexField
  field_type& fieldFiltered = 
    (0 == _vertexFilter) ? field : _vertexFilter->filter(field);
  field_type& fieldDimensioned = _dimension(fieldFiltered);
  
  _writer->writeVertexField(t, fieldDimensioned, mesh);
} // appendVertexField

// ----------------------------------------------------------------------
// Append finite-element cell field to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::OutputManager<mesh_type, field_type>::appendCellField(
				                     const PylithScalar t,
						     field_type& field,
						     const char* label,
						     const int labelId)
{ // appendCellField
  field_type& fieldFiltered = 
    (0 == _cellFilter) ? field : _cellFilter->filter(field, label, labelId);
  field_type& fieldDimensioned = _dimension(fieldFiltered);

  _writer->writeCellField(t, fieldDimensioned, label, labelId);
} // appendCellField

// ----------------------------------------------------------------------
// Dimension field.
template<typename mesh_type, typename field_type>
field_type&
pylith::meshio::OutputManager<mesh_type, field_type>::_dimension(field_type& fieldIn)
{ // _dimension
  if (1.0 == fieldIn.scale())
    return fieldIn;

  if (fieldIn.addDimensionOkay()) {
    fieldIn.dimensionalize();
    return fieldIn;
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
    
    if (0 == _fields)
      _fields = new topology::Fields<field_type>(fieldIn.mesh());
    
    if (!_fields->hasField(fieldName.c_str())) {
      ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
      logger.stagePush("Output");

      _fields->add(fieldName.c_str(), fieldIn.label());
      field_type& fieldOut = _fields->get(fieldName.c_str());
      fieldOut.cloneSection(fieldIn);
      fieldOut.vectorFieldType(fieldIn.vectorFieldType());
      fieldOut.scale(fieldIn.scale());

      logger.stagePop();
    } // if
    field_type& fieldOut = _fields->get(fieldName.c_str());
    fieldOut.copy(fieldIn);
    fieldOut.addDimensionOkay(true);
    fieldOut.dimensionalize();

    return fieldOut;
  } // if/else

  // Satisfy return value. Should never get this far.
  return fieldIn;
} // _dimension


// End of file 
