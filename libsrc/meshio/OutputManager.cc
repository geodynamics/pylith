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
// Set which vertex fields to output.
void
pylith::meshio::OutputManager::vertexFields(const char** names,
					    const int nfields)
{ // vertexFields
  _vertexFields.resize(nfields);
  for (int i=0; i < nfields; ++i)
    _vertexFields[i] = names[i];
} // vertexFields

// ----------------------------------------------------------------------
// Set which cell fields to output.
void
pylith::meshio::OutputManager::cellFields(const char** names,
					  const int nfields)
{ // cellFields
  _cellFields.resize(nfields);
  for (int i=0; i < nfields; ++i)
    _cellFields[i] = names[i];
} // cellFields

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

  const int nvfields = _vertexFields.size();
  for (int i=0; i < nvfields; ++i) {
    const ALE::Obj<real_section_type>& field = 
      fields->getReal(_vertexFields[i].c_str());
    
    // Create PETSc Vec for field values (if nec)
    // ADD STUFF HERE

    // Copy values from section to PETSc Vec
    // ADD STUFF HERE
    
    if (0 != _vertexFilter) {
      // Apply vertex filter
      // ADD STUFF HERE
    } // if
    //_writer->writeVertexField(t, fieldVec, _vertexFields[i], mesh);
  } // for

  const int ncfields = _cellFields.size();
  for (int i=0; i < ncfields; ++i) {
    const ALE::Obj<real_section_type>& field = 
      fields->getReal(_cellFields[i].c_str());
    
    // Create PETSc Vec for field values (if nec)
    // ADD STUFF HERE

    // Copy values from section to PETSc Vec
    // ADD STUFF HERE
    
    if (0 != _cellFilter) {
      // Apply vertex filter
      // ADD STUFF HERE
    } // if
    //_writer->writeCellField(t, fieldVec, _cellFields[i], mesh);
  } // for

  _writer->closeTimeStep();
} // writeFields


// End of file 
