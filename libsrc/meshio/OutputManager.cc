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

#include "DataWriter.hh" // USES DataWriter
#include "VertexFilter.hh" // USES VertexFilter
#include "CellFilter.hh" // USES CellFilter

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type>
pylith::meshio::OutputManager<mesh_type>::OutputManager(void) :
  _coordsys(0),
  _writer(0),
  _vertexFilter(0),
  _cellFilter(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type>
pylith::meshio::OutputManager<mesh_type>::~OutputManager(void)
{ // destructor
  _writer = 0; // :TODO: Use shared pointer
  _vertexFilter = 0; // :TODO: Use shared pointer
  _cellFilter = 0; // :TODO: Use shared pointer
  delete _coordsys; _coordsys = 0;
} // destructor  

// ----------------------------------------------------------------------
// Set coordinate system in output. The vertex fields in the output
template<typename mesh_type>
void
pylith::meshio::OutputManager<mesh_type>::coordsys(
				  const spatialdata::geocoords::CoordSys* cs)
{ // coordsys
  delete _coordsys; _coordsys = (0 != cs) ? cs->clone() : 0;
} // coordsys

// ----------------------------------------------------------------------
// Set writer to write data to file.
template<typename mesh_type>
void
pylith::meshio::OutputManager<mesh_type>::writer(
				       DataWriter<mesh_type>* const datawriter)
{ // writer
  _writer = datawriter; // :TODO: Use shared pointer
} // writer

// ----------------------------------------------------------------------
// Set filter for vertex data.
template<typename mesh_type>
void
pylith::meshio::OutputManager<mesh_type>::vertexFilter(
					VertexFilter<mesh_type>* const filter)
{ // vertexFilter
  _vertexFilter = filter; // :TODO: Use shared pointer
} // vertexFilter

// ----------------------------------------------------------------------
// Set filter for cell data.
template<typename mesh_type>
void
pylith::meshio::OutputManager<mesh_type>::cellFilter(
				         CellFilter<mesh_type>* const filter)
{ // cellFilter
  _cellFilter = filter; // :TODO: Use shared pointer
} // cellFilter

// ----------------------------------------------------------------------
// Prepare for output.
template<typename mesh_type>
void
pylith::meshio::OutputManager<mesh_type>::open(const mesh_type& mesh,
					       const int numTimeSteps,
					       const char* label,
					       const int labelId)
{ // open
  assert(0 != _writer);

  _writer->open(mesh, numTimeSteps, label, labelId);
} // open

// ----------------------------------------------------------------------
/// Close output files.
template<typename mesh_type>
void
pylith::meshio::OutputManager<mesh_type>::close(void)
{ // close
  assert(0 != _writer);
  _writer->close();
} // close

// ----------------------------------------------------------------------
// Setup file for writing fields at time step.
template<typename mesh_type>
void
pylith::meshio::OutputManager<mesh_type>::openTimeStep(const double t,
						       const mesh_type& mesh,
						       const char* label,
						       const int labelId)
{ // openTimeStep
  assert(0 != _writer);
  _writer->openTimeStep(t, mesh, label, labelId);
} // openTimeStep

// ----------------------------------------------------------------------
// End writing fields at time step.
template<typename mesh_type>
void
pylith::meshio::OutputManager<mesh_type>::closeTimeStep(void)
{ // closeTimeStep
  assert(0 != _writer);
  _writer->closeTimeStep();
} // closeTimeStep

// ----------------------------------------------------------------------
// Append finite-element vertex field to file.
template<typename mesh_type>
void
pylith::meshio::OutputManager<mesh_type>::appendVertexField(
			       const double t,
			       const topology::Field<mesh_type>& field)
{ // appendVertexField
  const topology::Field<mesh_type>& fieldFiltered = 
    (0 == _vertexFilter) ? field : _vertexFilter->filter(field);

  _writer->writeVertexField(t, fieldFiltered);
} // appendVertexField

// ----------------------------------------------------------------------
// Append finite-element cell field to file.
template<typename mesh_type>
void
pylith::meshio::OutputManager<mesh_type>::appendCellField(
				const double t,
				const topology::Field<mesh_type>& field,
				const char* label,
				const int labelId)
{ // appendCellField
  const topology::Field<mesh_type>& fieldFiltered = 
    (0 == _cellFilter) ? field : _cellFilter->filter(field, label, labelId);

  _writer->writeCellField(t, fieldFiltered, label, labelId);
} // appendCellField


// End of file 
