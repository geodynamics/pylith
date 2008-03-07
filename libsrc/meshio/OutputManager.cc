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
  delete _coordsys; _coordsys = 0;
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
				 const spatialdata::geocoords::CoordSys* csMesh,
				 const int numTimeSteps,
				 const char* label,
				 const int labelId)
{ // open
  assert(0 != _writer);

  _writer->open(mesh, csMesh, numTimeSteps, label, labelId);
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
// Setup file for writing fields at time step.
void
pylith::meshio::OutputManager::openTimeStep(
			     const double t,
			     const ALE::Obj<ALE::Mesh>& mesh,
			     const spatialdata::geocoords::CoordSys* csMesh,
			     const char* label,
			     const int labelId)
{ // openTimeStep
  assert(0 != _writer);
  _writer->openTimeStep(t, mesh, csMesh, label, labelId);
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
			       const VectorFieldEnum fieldType,
			       const ALE::Obj<ALE::Mesh>& mesh)
{ // appendVertexField
  assert(0 != name);

  VectorFieldEnum fieldTypeFiltered = fieldType;
  const ALE::Obj<real_section_type>& fieldFiltered = 
    (0 == _vertexFilter) ? 
    field : _vertexFilter->filter(&fieldTypeFiltered, field, mesh);

  _writer->writeVertexField(t, name, fieldFiltered, fieldTypeFiltered, mesh);
} // appendVertexField

// ----------------------------------------------------------------------
// Append finite-element cell field to file.
void
pylith::meshio::OutputManager::appendCellField(
				const double t,
				const char* name,
				const ALE::Obj<real_section_type>& field,
				const VectorFieldEnum fieldType,
				const ALE::Obj<ALE::Mesh>& mesh,
				const char* label,
				const int labelId)
{ // appendCellField
  assert(0 != name);

  VectorFieldEnum fieldTypeFiltered = fieldType;
  const ALE::Obj<real_section_type>& fieldFiltered = 
    (0 == _cellFilter) ? 
    field : _cellFilter->filter(&fieldTypeFiltered, field, 
				mesh, label, labelId);

  _writer->writeCellField(t, name, fieldFiltered, fieldTypeFiltered,
			  mesh, label, labelId);
} // appendCellField


// End of file 
