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

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::DataWriterVTK::DataWriterVTK(void) :
  _timeConstant(1.0),
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
  _timeConstant(w._timeConstant),
  _filename(w._filename),
  _timeFormat(w._timeFormat),
  _viewer(0),
  _wroteVertexHeader(w._wroteVertexHeader),
  _wroteCellHeader(w._wroteCellHeader)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Set value used to normalize time stamp in name of VTK file.
void
pylith::meshio::DataWriterVTK::timeConstant(const double value)
{ // timeConstant
  if (value <= 0.0) {
    std::ostringstream msg;
    msg << "Time used to normalize time stamp in VTK data files must be "
	<< "positive.\nCurrent value is " << value << ".";
    throw std::runtime_error(msg.str());
  } // if
  _timeConstant = value;
} // timeConstant

// ----------------------------------------------------------------------
// Prepare file for data at a new time step.
void
pylith::meshio::DataWriterVTK::openTimeStep(
			       const double t,
			       const ALE::Obj<Mesh>& mesh,
			       const spatialdata::geocoords::CoordSys* csMesh,
			       const char* label,
			       const int labelId)
{ // openTimeStep
  assert(!mesh.isNull());
  assert(0 != csMesh);

  try {
    PetscErrorCode err;

    const std::string& filename = _vtkFilename(t);

    err = PetscViewerCreate(mesh->comm(), &_viewer);
    err = PetscViewerSetType(_viewer, PETSC_VIEWER_ASCII);
    err = PetscViewerSetFormat(_viewer, PETSC_VIEWER_ASCII_VTK);
    err = PetscViewerFileSetName(_viewer, filename.c_str());
    if (err)
      throw std::runtime_error("Could not open VTK file.");
    
    err = VTKViewer::writeHeader(_viewer);
    //std::cout << "Wrote header for " << filename << std::endl;
    err = VTKViewer::writeVertices(mesh, _viewer);
    //std::cout << "Wrote vertices for " << filename << std::endl;
    if (0 == label)
      err = VTKViewer::writeElements(mesh, _viewer);
    else {
      const std::string labelName = 
	(mesh->hasLabel("censored depth")) ? "censored depth" : "depth";
      err = VTKViewer::writeElements(mesh, label, labelId, labelName, 0, _viewer);      
    } // if
    if (err)
      throw std::runtime_error("Could not write topology.");
    //std::cout << "Wrote elements for " << filename << std::endl;

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
				       const ALE::Obj<Mesh>& mesh)
{ // writeVertexField
  assert(0 != name);
  assert(!mesh.isNull());
  assert(!field.isNull());

  try {
    const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
    assert(!vertices.isNull());
    int rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    const std::string labelName = 
      (mesh->hasLabel("censored depth")) ? "censored depth" : "depth";
    const ALE::Obj<Mesh::numbering_type>& numbering =
      mesh->getFactory()->getNumbering(mesh, labelName, 0);
    assert(!numbering.isNull());

    const int localFiberDim = 
      field->getFiberDimension(*mesh->getLabelStratum(labelName, 0)->begin());
    int fiberDim;
    MPI_Allreduce((void *) &localFiberDim, (void *) &fiberDim, 1, MPI_INT, MPI_MAX,
		  mesh->comm());
    assert(fiberDim > 0);
    const int enforceDim = (fieldType != VECTOR_FIELD) ? fiberDim : 3;

    PetscErrorCode err = 0;

    if (!_wroteVertexHeader) {
      err = PetscViewerASCIIPrintf(_viewer, "POINT_DATA %d\n", 
						  numbering->getGlobalSize());
      if (err)
	throw std::runtime_error("Could not write VTK point data header.");
      _wroteVertexHeader = true;
    } // if

    err = VTKViewer::writeField(field, name, fiberDim, numbering, _viewer, 
				enforceDim);
    if (err)
      throw std::runtime_error("Coult not write vertex field.");

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
				       const ALE::Obj<Mesh>& mesh,
				       const char* label,
				       const int labelId)
{ // writeCellField
  assert(0 != name);
  assert(!mesh.isNull());
  assert(!field.isNull());

  try {
    const ALE::Obj<Mesh::label_sequence>& cells = (0 == label) ?
      mesh->heightStratum(0) :
      mesh->getLabelStratum(label, labelId);
    assert(!cells.isNull());
    int rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Correctly handle boundary and fault meshes
    //   Cannot just use mesh->depth() because boundaries report the wrong thing
    const int cellDepth = (mesh->depth() == -1) ? -1 : 1;
    const int depth = (0 == label) ? cellDepth : labelId;
    const std::string labelName = (0 == label) ?
      ((mesh->hasLabel("censored depth")) ? "censored depth" : "depth") : label;
    const ALE::Obj<Mesh::numbering_type>& numbering = 
      mesh->getFactory()->getNumbering(mesh, labelName, depth);
    assert(!numbering.isNull());
    const int localFiberDim = 
      field->getFiberDimension(*mesh->getLabelStratum(labelName, depth)->begin());
    int fiberDim;
    MPI_Allreduce((void *) &localFiberDim, (void *) &fiberDim, 1, MPI_INT, MPI_MAX, 
		  mesh->comm());
    assert(fiberDim > 0);
    const int enforceDim = (fieldType != VECTOR_FIELD) ? fiberDim : 3;

    PetscErrorCode err = 0;

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

// ----------------------------------------------------------------------
// Generate filename for VTK file.
std::string
pylith::meshio::DataWriterVTK::_vtkFilename(const double t) const
{ // _vtkFilename
  std::ostringstream filename;
  const int indexExt = _filename.find(".vtk");
  if (_numTimeSteps > 0) {
    // If data with multiple time steps, then add time stamp to filename
    char sbuffer[256];
    sprintf(sbuffer, _timeFormat.c_str(), t/_timeConstant);
    std::string timestamp(sbuffer);
    const int pos = timestamp.find(".");
    if (pos >0 && pos != timestamp.length())
      timestamp.erase(pos, 1);
    filename
      << std::string(_filename, 0, indexExt) << "_t" << timestamp << ".vtk";
  } else
    filename
      << std::string(_filename, 0, indexExt) << "_info.vtk";

  return std::string(filename.str());
} // _vtkFilename


// End of file 
