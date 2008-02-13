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
  _viewer(0),
  _wroteVertexHeader(false),
  _wroteCellHeader(false)
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
  _viewer(0),
  _wroteVertexHeader(w._wroteVertexHeader),
  _wroteCellHeader(w._wroteCellHeader)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Prepare file for data at a new time step.
void
pylith::meshio::DataWriterVTK::openTimeStep(
			       const double t,
			       const ALE::Obj<ALE::Mesh>& mesh,
			       const spatialdata::geocoords::CoordSys* csMesh,
			       const char* label,
			       const int labelId)
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
    if (0 == label)
      err = VTKViewer::writeElements(mesh, _viewer);
    else
      err = VTKViewer::writeElements(mesh, label, labelId, "censored depth", 0, _viewer);      
    if (err)
      throw std::runtime_error("Could not write topology.");

    _wroteVertexHeader = false;
    _wroteCellHeader = false;
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
  _wroteVertexHeader = false;
  _wroteCellHeader = false;
} // closeTimeStep

// ----------------------------------------------------------------------
// Write field over vertices to file.
void
pylith::meshio::DataWriterVTK::writeVertexField(
				       const double t,
				       const char* name,
				       const ALE::Obj<real_section_type>& field,
				       const VectorFieldEnum fieldType,
				       const ALE::Obj<ALE::Mesh>& mesh)
{ // writeVertexField
  assert(0 != name);

  try {
    const std::string labelName = 
      (mesh->hasLabel("censored depth")) ? "censored depth" : "depth";
    const ALE::Obj<Mesh::numbering_type>& numbering =
      mesh->getFactory()->getNumbering(mesh, labelName, 0);

    const int fiberDim = 
      field->getFiberDimension(*mesh->getLabelStratum(labelName, 0)->begin());
    const int enforceDim = (fieldType != VECTOR_FIELD) ? fiberDim : 3;

    if (!_wroteVertexHeader) {
      PetscErrorCode err = PetscViewerASCIIPrintf(_viewer, "POINT_DATA %d\n", 
						  numbering->getGlobalSize());
      if (err)
	throw std::runtime_error("Could not write VTK point data header.");
      _wroteVertexHeader = true;
    } // if

    VTKViewer::writeField(field, name, fiberDim, numbering, _viewer, 
			  enforceDim);

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
				       const VectorFieldEnum fieldType,
				       const ALE::Obj<ALE::Mesh>& mesh)
{ // writeCellField
  assert(0 != name);

  try {
    std::string    labelName;
    PetscErrorCode err = 0;

    if (mesh->hasLabel("censored depth")) {
      labelName = "censored depth";
    } else {
      labelName = "depth";
    }

    const ALE::Obj<Mesh::numbering_type>& numbering = 
      mesh->getFactory()->getNumbering(mesh, labelName, mesh->depth());
    const int fiberDim = 
      field->getFiberDimension(*mesh->getLabelStratum(labelName, 
						      mesh->depth())->begin());
    const int enforceDim = (fieldType != VECTOR_FIELD) ? fiberDim : 3;

    if (!_wroteCellHeader) {
      err = PetscViewerASCIIPrintf(_viewer, "CELL_DATA %d\n", 
				   numbering->getGlobalSize());
      if (err)
	throw std::runtime_error("Could not write VTK point data header.");
      _wroteCellHeader = true;
    } // if

    VTKViewer::writeField(field, name, fiberDim, numbering, _viewer, 
			  enforceDim);

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
