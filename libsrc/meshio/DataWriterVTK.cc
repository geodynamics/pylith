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

#include "DataWriterVTK.hh" // implementation of class methods

#include <petscmesh_viewers.hh> // USES VTKViewer

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::DataWriterVTK::DataWriterVTK(void) :
  _filename("output.vtk"),
  _timeFormat("%f"),
  _viewer(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriterVTK::~DataWriterVTK(void)
{ // destructor
  if (0 != _viewer)
    PetscViewerDestroy(_viewer);
  _viewer = 0;
} // destructor  

// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::DataWriterVTK::DataWriterVTK(const DataWriterVTK& w) :
  DataWriter(w),
  _filename(w._filename),
  _timeFormat(w._timeFormat),
  _viewer(0)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Prepare file for data at a new time step.
void
pylith::meshio::DataWriterVTK::openTimeStep(
			       const double t,
			       const ALE::Obj<ALE::Mesh>& mesh,
			       const spatialdata::geocoords::CoordSys* csMesh)
{ // openTimeStep
  assert(!mesh.isNull());
  assert(0 != csMesh);

  try {
    PetscErrorCode err;

    std::ostringstream buffer;
    const int indexExt = _filename.find(".vtk");
    if (_numTimeSteps > 0) {
      // If data with multiple time steps, then add time stamp to filename
      char sbuffer[256];
      sprintf(sbuffer, _timeFormat.c_str(), t);
      std::string timestamp(sbuffer);
      const int pos = timestamp.find(".");
      if (pos != timestamp.length())
	timestamp.erase(pos, 1);
      buffer
	<< std::string(_filename, 0, indexExt) << "_t" << timestamp << ".vtk";
    } else
      buffer
	<< std::string(_filename, 0, indexExt) << "_info.vtk";

    err = PetscViewerCreate(mesh->comm(), &_viewer);
    err = PetscViewerSetType(_viewer, PETSC_VIEWER_ASCII);
    err = PetscViewerSetFormat(_viewer, PETSC_VIEWER_ASCII_VTK);
    err = PetscViewerFileSetName(_viewer, buffer.str().c_str());
    if (err)
      throw std::runtime_error("Could not open VTK file.");
    
    err = VTKViewer::writeHeader(_viewer);
    err = VTKViewer::writeVertices(mesh, _viewer);
    err = VTKViewer::writeElements(mesh, _viewer);
    if (err)
      throw std::runtime_error("Could not write topology.");
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while preparing for writing data to VTK file "
	<< _filename << " at time " << t << ".\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (const ALE::Exception& err) {
    std::ostringstream msg;
    msg << "Error while preparing for writing data to VTK file "
	<< _filename << " at time " << t << ".\n" << err.msg();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Unknown error while preparing for writing data to VTK file "
	<< _filename << " at time " << t << ".\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // openTimeStep

// ----------------------------------------------------------------------
/// Cleanup after writing data for a time step.
void
pylith::meshio::DataWriterVTK::closeTimeStep(void)
{ // closeTimeStep
  PetscViewerDestroy(_viewer); _viewer = 0;
} // closeTimeStep

// ----------------------------------------------------------------------
// Write field over vertices to file.
void
pylith::meshio::DataWriterVTK::writeVertexField(
				       const double t,
				       const char* name,
				       const ALE::Obj<real_section_type>& field,
				       const FieldEnum fieldType,
				       const ALE::Obj<ALE::Mesh>& mesh)
{ // writeVertexField
  assert(0 != name);

  try {
    std::ostringstream buffer;
    buffer.str("");
    char timestamp[256];
    sprintf(timestamp, _timeFormat.c_str(), t);
    buffer << name << "_t" << timestamp;

    const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
    const int fiberDim = (fieldType != VECTOR_FIELD) ? 
      field->getFiberDimension(*vertices->begin()) : 3;

    PetscErrorCode err = SectionView_Sieve_Ascii(mesh, field, 
						 buffer.str().c_str(), 
						 _viewer, fiberDim);
    if (err)
      throw std::runtime_error("Could not write vertex data.");
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << name << "' at time " 
	<< t << " to VTK file '" << _filename << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << name << "' at time " 
	<< t << " to VTK file '" << _filename << "'.\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // writeVertexField

// ----------------------------------------------------------------------
// Write field over cells to file.
void
pylith::meshio::DataWriterVTK::writeCellField(
				       const double t,
				       const char* name,
				       const ALE::Obj<real_section_type>& field,
				       const FieldEnum fieldType,
				       const ALE::Obj<ALE::Mesh>& mesh)
{ // writeCellField
  assert(0 != name);

  try {
    PetscErrorCode err = 0;

    std::ostringstream buffer;
    buffer.str("");
    char timestamp[256];
    sprintf(timestamp, _timeFormat.c_str(), t);
    buffer << name << "_t" << timestamp;

    err = PetscViewerPushFormat(_viewer, PETSC_VIEWER_ASCII_VTK_CELL);
    
    // Get fiber dimension of first cell
    const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);
    const int fiberDim = (fieldType != VECTOR_FIELD) ? 
      field->getFiberDimension(*cells->begin()) : 3;
    err = SectionView_Sieve_Ascii(mesh, field, buffer.str().c_str(), 
				  _viewer, fiberDim);
    if (err)
      throw std::runtime_error("Could not write cell data.");   
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << name << "' at time " 
	<< t << " to VTK file '" << _filename << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << name << "' at time " 
	<< t << " to VTK file '" << _filename << "'.\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // writeCellField


// End of file 
