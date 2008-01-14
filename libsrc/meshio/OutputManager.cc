// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "OutputManager.hh" // implementation of class methods

#include "DataWriter.hh" // USES DataWriter
#include "VertexFilter.hh" // USES VertexFilter
#include "CellFilter.hh" // USES CellFilter
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputManager::OutputManager(void) :
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
  delete _writer; _writer = 0;
  delete _vertexFilter; _vertexFilter = 0;
  delete _cellFilter; _cellFilter = 0;
} // destructor  

// ----------------------------------------------------------------------
// Set coordinate system in output. The vertex fields in the output
void
pylith::meshio::OutputManager::coordsys(const spatialdata::geocoords::CoordSys* cs)
{ // coordsys
  delete _coordsys; _coordsys = (0 != cs) ? cs->clone() : 0;
} // coordsys

// ----------------------------------------------------------------------
// Set writer to write data to file.
void
pylith::meshio::OutputManager::writer(const DataWriter* datawriter)
{ // writer
  delete _writer; _writer = (0 != datawriter) ? datawriter->clone() : 0;
} // writer

// ----------------------------------------------------------------------
// Set which vertex fields to output.
void
pylith::meshio::OutputManager::vertexFields(const char** names,
					    const char** labels,
					    const int numFields)
{ // vertexFields
  assert( (0 == numFields && 0 == names && 0 == labels) ||
	  (0 < numFields && 0 != names && 0 != labels) );

  _vertexFields.clear();
  for (int iField=0; iField < numFields; ++iField)
    _vertexFields[names[iField]] = _vertexFields[labels[iField]];
} // vertexFields

// ----------------------------------------------------------------------
// Set which cell fields to output.
void
pylith::meshio::OutputManager::cellFields(const char** names,
					  const char** labels,
					  const int numFields)
{ // cellFields
  assert( (0 == numFields && 0 == names && 0 == labels) ||
	  (0 < numFields && 0 != names && 0 != labels) );

  _cellFields.clear();
  for (int iField=0; iField < numFields; ++iField)
    _cellFields[names[iField]] = _cellFields[labels[iField]];
} // cellFields

// ----------------------------------------------------------------------
// Get vertex fields to output.
const pylith::meshio::OutputManager::map_names_type&
pylith::meshio::OutputManager::vertexFields(void) const
{ // vertexFields
  return _vertexFields;
} // vertexFields

// ----------------------------------------------------------------------
// Get cell fields to output.
const pylith::meshio::OutputManager::map_names_type&
pylith::meshio::OutputManager::cellFields(void) const
{ // cellFields
  return _cellFields;
} // cellFields

// ----------------------------------------------------------------------
// Set filter for vertex data.
void
pylith::meshio::OutputManager::vertexFilter(const VertexFilter* filter)
{ // vertexFilter
  delete _vertexFilter; _vertexFilter = (0 != filter) ? filter->clone() : 0;
} // vertexFilter

// ----------------------------------------------------------------------
// Set filter for cell data.
void
pylith::meshio::OutputManager::cellFilter(const CellFilter* filter)
{ // cellFilter
  delete _cellFilter; _cellFilter = (0 != filter) ? filter->clone() : 0;
} // cellFilter

// ----------------------------------------------------------------------
// Prepare for output.
void
pylith::meshio::OutputManager::open(
				 const ALE::Obj<ALE::Mesh>& mesh,
				 const spatialdata::geocoords::CoordSys* csMesh)
{ // open
  assert(0 != _writer);

  _writer->open(mesh, csMesh);
} // open

// ----------------------------------------------------------------------
/// Close output files.
void
pylith::meshio::OutputManager::close(void)
{ // close
  assert(0 != _writer);
  _writer->close();
} // close

// ----------------------------------------------------------------------
// Write finite-element fields to file.
void
pylith::meshio::OutputManager::writeFields(
				const double t,
				topology::FieldsManager* const fields,
				const ALE::Obj<ALE::Mesh>& mesh,
				const spatialdata::geocoords::CoordSys* csMesh)
{ // writeFields
  assert(0 != _writer);
  assert(0 != fields);

  _writer->openTimeStep(t, mesh, csMesh);
  
  for (map_names_type::iterator f_iter=_vertexFields.begin();
       f_iter != _vertexFields.end();
       ++f_iter) {
    const char* fieldName = f_iter->first.c_str();
    const char* fieldLabel = f_iter->second.c_str();
    const ALE::Obj<real_section_type>& field = fields->getReal(fieldLabel);

    const ALE::Obj<real_section_type>& fieldFiltered = 
      (0 != _vertexFilter) ? field : _vertexFilter->filter(field, mesh);

    _writer->writeVertexField(t, fieldName, fieldFiltered, mesh);
  } // for

  for (map_names_type::iterator f_iter=_cellFields.begin();
       f_iter != _cellFields.end();
       ++f_iter) {
    const char* fieldName = f_iter->first.c_str();
    const char* fieldLabel = f_iter->second.c_str();
    const ALE::Obj<real_section_type>& field = fields->getReal(fieldLabel);
    
    const ALE::Obj<real_section_type>& fieldFiltered = 
      (0 != _cellFilter) ? field : _cellFilter->filter(field, mesh);

    _writer->writeCellField(t, fieldName, fieldFiltered, mesh);
  } // for

  _writer->closeTimeStep();
} // writeFields

// ----------------------------------------------------------------------
// Setup file for writing fields at time step.
void
pylith::meshio::OutputManager::openTimeStep(const double t,
	     const ALE::Obj<ALE::Mesh>& mesh,
	     const spatialdata::geocoords::CoordSys* csMesh)
{ // openTimeStep
  assert(0 != _writer);
  _writer->openTimeStep(t, mesh, csMesh);
} // openTimeStep

// ----------------------------------------------------------------------
// End writing fields at time step.
void
pylith::meshio::OutputManager::closeTimeStep(void)
{ // closeTimeStep
  assert(0 != _writer);
  _writer->closeTimeStep();
} // closeTimeStep

// ----------------------------------------------------------------------
// Append finite-element vertex field to file.
void
pylith::meshio::OutputManager::appendVertexField(
			       const double t,
			       const char* name,
			       const ALE::Obj<real_section_type>& field,
			       const ALE::Obj<ALE::Mesh>& mesh)
{ // appendVertexField
  assert(0 != name);

  const ALE::Obj<real_section_type>& fieldFiltered = 
    (0 != _vertexFilter) ? field : _vertexFilter->filter(field, mesh);

  _writer->writeVertexField(t, name, fieldFiltered, mesh);
} // appendVertexField

// ----------------------------------------------------------------------
// Append finite-element cell field to file.
void
pylith::meshio::OutputManager::appendCellField(
				const double t,
				const char* name,
				const ALE::Obj<real_section_type>& field,
				const ALE::Obj<ALE::Mesh>& mesh)
{ // appendCellField
  assert(0 != name);

  const ALE::Obj<real_section_type>& fieldFiltered = 
    (0 != _cellFilter) ? field : _cellFilter->filter(field, mesh);

  _writer->writeCellField(t, name, fieldFiltered, mesh);
} // appendCellField


// End of file 
