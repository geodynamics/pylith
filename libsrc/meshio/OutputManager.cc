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

#include "DataWriter.hh" // HOLDSA DataWriter
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
#if 0 // TEMPORARY
#include "VertexFilter.hh" // HOLDS VertexFilter
#include "CellFilter.hh" // HOLDS CellFilter
#endif // TEMPORARY

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputManager::OutputManager(void) :
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
#if 0 // TEMPORARY
  delete _vertexFilter; _vertexFilter = 0;
  delete _cellFilter; _cellFilter = 0;
#endif // TEMPORARY
} // destructor  

// ----------------------------------------------------------------------
// Set filter for vertex data.
void
pylith::meshio::OutputManager::vertexFilter(const VertexFilter* filter)
{ // vertexFilter
#if 0 // TEMPORARY
  delete _vertexFilter; _vertexFilter = (0 != filter) ? filter->clone() : 0;
#endif // TEMPORARY
} // vertexFilter

// ----------------------------------------------------------------------
// Set filter for cell data.
void
pylith::meshio::OutputManager::cellFilter(const CellFilter* filter)
{ // cellFilter
#if 0 // TEMPORARY
  delete _cellFilter; _cellFilter = (0 != filter) ? filter->clone() : 0;
#endif // TEMPORARY
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

#if 0
    const ALE::Obj<real_section_type>& fieldFiltered = 
      (0 != _vertexFilter) ? field : _vertexFilter->filter(field);

    _writer->writeVertexField(t, fieldFiltered, fieldName, mesh);
#endif
  } // for

  for (map_names_type::iterator f_iter=_cellFields.begin();
       f_iter != _cellFields.end();
       ++f_iter) {
    const char* fieldName = f_iter->first.c_str();
    const char* fieldLabel = f_iter->second.c_str();
    const ALE::Obj<real_section_type>& field = fields->getReal(fieldLabel);
    
#if 0
    const ALE::Obj<real_section_type>& fieldFiltered = 
      (0 != _cellFilter) ? field : _cellFilter->filter(field);

    _writer->writeCellField(t, fieldFiltered, fieldName, mesh);
#endif
  } // for

  _writer->closeTimeStep();
} // writeFields


// End of file 
