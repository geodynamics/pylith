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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include <petscdmmesh_viewers.hh> // USES VTKViewer

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterVTK<mesh_type,field_type>::DataWriterVTK(void) :
  _timeConstant(1.0),
  _filename("output.vtk"),
  _timeFormat("%f"),
  _viewer(0),
  _wroteVertexHeader(false),
  _wroteCellHeader(false),
  _precision(6)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterVTK<mesh_type,field_type>::~DataWriterVTK(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Set precision of floating point values in output.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type,field_type>::precision(const int value)
{ // precision
  if (value <= 0) {
    std::ostringstream msg;
    msg << "Floating point precision (" << value << ") must be positive.";
    throw std::runtime_error(msg.str());
  } // if

  _precision = value;
} // precision

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type, field_type>::deallocate(void)
{ // deallocate
  DataWriter<mesh_type, field_type>::deallocate();

  PetscErrorCode err = 0;
  err = PetscViewerDestroy(&_viewer);CHECK_PETSC_ERROR(err);
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor.
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterVTK<mesh_type,field_type>::DataWriterVTK(const DataWriterVTK<mesh_type, field_type>& w) :
  DataWriter<mesh_type, field_type>(w),
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
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type,field_type>::timeConstant(const PylithScalar value)
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
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type,field_type>::openTimeStep(const PylithScalar t,
						       const mesh_type& mesh,
						       const char* label,
						       const int labelId)
{ // openTimeStep

  try {
    PetscErrorCode err = 0;
    
    const std::string& filename = _vtkFilename(t);

    err = PetscViewerCreate(mesh.comm(), &_viewer);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerSetType(_viewer, PETSCVIEWERASCII);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerSetFormat(_viewer, PETSC_VIEWER_ASCII_VTK);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerFileSetName(_viewer, filename.c_str());
    CHECK_PETSC_ERROR(err);

    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();

    err = VTKViewer::writeHeader(sieveMesh, _viewer);
    CHECK_PETSC_ERROR(err);
    //std::cout << "Wrote header for " << filename << std::endl;
    err = VTKViewer::writeVertices(sieveMesh, _viewer);
    CHECK_PETSC_ERROR(err);
    //std::cout << "Wrote vertices for " << filename << std::endl;

    const int cellDepth = (sieveMesh->depth() == -1) ? -1 : 1;
    const int depth = (!label) ? cellDepth : labelId;
    const std::string cLabelName = (!label) ?
      ((sieveMesh->hasLabel("censored depth")) ?
       "censored depth" : "depth") : label;
    const std::string vLabelName = 
      (sieveMesh->hasLabel("censored depth")) ? "censored depth" : "depth";
    err = VTKViewer::writeElements(sieveMesh, cLabelName, depth, vLabelName, 0,
				   _viewer); CHECK_PETSC_ERROR(err);
    //std::cout << "Wrote elements for " << filename << std::endl;

    _wroteVertexHeader = false;
    _wroteCellHeader = false;
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while preparing for writing data to VTK file "
	<< _vtkFilename(t) << " at time " << t << ".\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Unknown error while preparing for writing data to VTK file "
	<< _vtkFilename(t) << " at time " << t << ".\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // openTimeStep

// ----------------------------------------------------------------------
/// Cleanup after writing data for a time step.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type,field_type>::closeTimeStep(void)
{ // closeTimeStep
  PetscViewerDestroy(&_viewer); _viewer = 0;
  _wroteVertexHeader = false;
  _wroteCellHeader = false;
} // closeTimeStep

// ----------------------------------------------------------------------
// Write field over vertices to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type,field_type>::writeVertexField(
				            const PylithScalar t,
					    field_type& field,
					    const mesh_type& mesh)
{ // writeVertexField
  typedef typename mesh_type::SieveMesh SieveMesh;
  typedef typename field_type::Mesh::RealSection RealSection;

  try {
    const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
    assert(!sieveMesh.isNull());

    const std::string labelName = 
      (sieveMesh->hasLabel("censored depth")) ? "censored depth" : "depth";
    const ALE::Obj<typename SieveMesh::numbering_type>& numbering =
      sieveMesh->getFactory()->getNumbering(sieveMesh, labelName, 0);
    assert(!numbering.isNull());

    const ALE::Obj<RealSection>& section = field.section();
    assert(!section.isNull());
    assert(!sieveMesh->getLabelStratum(labelName, 0).isNull());
    
    int fiberDimLocal = 
      (sieveMesh->getLabelStratum(labelName, 0)->size() > 0) ? 
      section->getFiberDimension(*sieveMesh->getLabelStratum(labelName, 0)->begin()) : 0;
    int fiberDim = 0;
    MPI_Allreduce(&fiberDimLocal, &fiberDim, 1, MPI_INT, MPI_MAX,
		  field.mesh().comm());
    assert(fiberDim > 0);
    const int enforceDim =
      (field.vectorFieldType() != topology::FieldBase::VECTOR) ? fiberDim : 3;

    PetscErrorCode err = 0;
    if (!_wroteVertexHeader) {
      err = PetscViewerASCIIPrintf(_viewer, "POINT_DATA %d\n", 
				   numbering->getGlobalSize());
      CHECK_PETSC_ERROR(err);
      _wroteVertexHeader = true;
    } // if

    err = VTKViewer::writeField(section, field.label(), fiberDim, numbering,
				_viewer, enforceDim, _precision);
    CHECK_PETSC_ERROR(err);
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to VTK file '" << _vtkFilename(t) << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to VTK file '" << _vtkFilename(t) << "'.\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // writeVertexField

// ----------------------------------------------------------------------
// Write field over cells to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type,field_type>::writeCellField(
				       const PylithScalar t,
				       field_type& field,
				       const char* label,
				       const int labelId)
{ // writeCellField
  typedef typename field_type::Mesh::SieveMesh SieveMesh;
  typedef typename field_type::Mesh::RealSection RealSection;

  try {
    PetscErrorCode err = 0;

    // Correctly handle boundary and fault meshes
    //   Cannot just use mesh->depth() because boundaries report the wrong thing
    const ALE::Obj<SieveMesh>& sieveMesh = field.mesh().sieveMesh();
    assert(!sieveMesh.isNull());
    int cellDepthLocal = (sieveMesh->depth() == -1) ? -1 : 1;
    int cellDepth = 0;
    err = MPI_Allreduce(&cellDepthLocal, &cellDepth, 1, MPI_INT, MPI_MAX, 
			sieveMesh->comm());CHECK_PETSC_ERROR(err);
    const int depth = (!label) ? cellDepth : labelId;
    const std::string labelName = (!label) ?
      ((sieveMesh->hasLabel("censored depth")) ?
       "censored depth" : "depth") : label;
    assert(!sieveMesh->getFactory().isNull());
    const ALE::Obj<typename SieveMesh::numbering_type>& numbering = 
      sieveMesh->getFactory()->getNumbering(sieveMesh, labelName, depth);
    assert(!numbering.isNull());
    const ALE::Obj<RealSection>& section = field.section();
    assert(!section.isNull());

    assert(!sieveMesh->getLabelStratum(labelName, depth).isNull());
    int fiberDimLocal = 
      (sieveMesh->getLabelStratum(labelName, depth)->size() > 0) ? 
      section->getFiberDimension(*sieveMesh->getLabelStratum(labelName, depth)->begin()) : 0;
    int fiberDim = 0;
    MPI_Allreduce(&fiberDimLocal, &fiberDim, 1, MPI_INT, MPI_MAX,
		  field.mesh().comm());
    assert(fiberDim > 0);

    const int enforceDim =
      (field.vectorFieldType() != topology::FieldBase::VECTOR) ? fiberDim : 3;

    if (!_wroteCellHeader) {
      err = PetscViewerASCIIPrintf(_viewer, "CELL_DATA %d\n", 
				   numbering->getGlobalSize());
      CHECK_PETSC_ERROR(err);
      _wroteCellHeader = true;
    } // if

    VTKViewer::writeField(section, field.label(), fiberDim, numbering,
			  _viewer, enforceDim, _precision);
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to VTK file '" << _vtkFilename(t) << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to VTK file '" << _vtkFilename(t) << "'.\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // writeCellField

// ----------------------------------------------------------------------
// Generate filename for VTK file.
template<typename mesh_type, typename field_type>
std::string
pylith::meshio::DataWriterVTK<mesh_type,field_type>::_vtkFilename(const PylithScalar t) const
{ // _vtkFilename
  std::ostringstream filename;
  const int indexExt = _filename.find(".vtk");
  const int numTimeSteps = DataWriter<mesh_type, field_type>::_numTimeSteps;
  if (numTimeSteps > 0) {
    // If data with multiple time steps, then add time stamp to filename
    char sbuffer[256];
    sprintf(sbuffer, _timeFormat.c_str(), t/_timeConstant);
    std::string timestamp(sbuffer);
    const int pos = timestamp.find(".");
    if (pos > 0 && pos != timestamp.length())
      timestamp.erase(pos, 1);
    filename
      << std::string(_filename, 0, indexExt) << "_t" << timestamp << ".vtk";
  } else
    filename
      << std::string(_filename, 0, indexExt) << "_info.vtk";

  return std::string(filename.str());
} // _vtkFilename


// End of file 
