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

#include "pylith/topology/Stratum.hh" // USES StratumIS

#include <petscdmmesh_viewers.hh> // USES VTKViewer
#include <petscdmplex.h>

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

extern PetscErrorCode DMPlexVTKWriteAll(PetscObject odm, PetscViewer viewer);

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterVTK<mesh_type,field_type>::DataWriterVTK(void) :
  _timeConstant(1.0),
  _filename("output.vtk"),
  _timeFormat("%f"),
  _viewer(NULL),
  _dm(NULL),
  _isOpen(false),
  _isOpenTimeStep(false),
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
// Deallocate PETSc and local data structures.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type, field_type>::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  closeTimeStep(); // Insure time step is closed.
  close(); // Insure clean up.
  DataWriter<mesh_type, field_type>::deallocate();

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor.
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterVTK<mesh_type,field_type>::DataWriterVTK(const DataWriterVTK<mesh_type, field_type>& w) :
  DataWriter<mesh_type, field_type>(w),
  _timeConstant(w._timeConstant),
  _filename(w._filename),
  _timeFormat(w._timeFormat),
  _viewer(NULL),
  _dm(NULL),
  _isOpen(w._isOpen),
  _isOpenTimeStep(w._isOpenTimeStep),
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
  PYLITH_METHOD_BEGIN;

  if (value <= 0.0) {
    std::ostringstream msg;
    msg << "Time used to normalize time stamp in VTK data files must be "
	<< "positive.\nCurrent value is " << value << ".";
    throw std::runtime_error(msg.str());
  } // if
  _timeConstant = value;

  PYLITH_METHOD_END;
} // timeConstant

// ----------------------------------------------------------------------
// Set precision of floating point values in output.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type,field_type>::precision(const int value)
{ // precision
  PYLITH_METHOD_BEGIN;

  if (value <= 0) {
    std::ostringstream msg;
    msg << "Floating point precision (" << value << ") must be positive.";
    throw std::runtime_error(msg.str());
  } // if

  _precision = value;

  PYLITH_METHOD_END;
} // precision

// ----------------------------------------------------------------------
// Prepare for writing files.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type,field_type>::open(const mesh_type& mesh,
							  const int numTimeSteps,
							  const char* label,
							  const int labelId)
{ // open
  PYLITH_METHOD_BEGIN;

  DataWriter<mesh_type, field_type>::open(mesh, numTimeSteps, label, labelId);

  // Save handle for actions required in closeTimeStep() and close();
  PetscErrorCode err = 0;
  err = DMDestroy(&_dm);PYLITH_CHECK_ERROR(err);
  _dm = mesh.dmMesh();assert(_dm);
  err = PetscObjectReference((PetscObject) _dm);PYLITH_CHECK_ERROR(err);

  // Create VTK label in DM: Cleared in close().
  if (label) {
    topology::StratumIS cellsIS(_dm, label, labelId);
    const PetscInt ncells = cellsIS.size();
    const PetscInt* cells = cellsIS.points();

    for (PetscInt c=0; c < ncells; ++c) {
      err = DMPlexSetLabelValue(_dm, "vtk", cells[c], 1);PYLITH_CHECK_ERROR(err);
    } // for

  } // if

  _isOpen = true;

  PYLITH_METHOD_END;
} // open

// ----------------------------------------------------------------------
// Close output files.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type,field_type>::close(void)
{ // close
  PYLITH_METHOD_BEGIN;

  if (_isOpen) {
    assert(_dm);
    PetscBool hasLabel = PETSC_FALSE;
    PetscErrorCode err = DMPlexHasLabel(_dm, "vtk", &hasLabel);PYLITH_CHECK_ERROR(err);
    if (hasLabel) {
      err = DMPlexClearLabelStratum(_dm, "vtk", 1);PYLITH_CHECK_ERROR(err);
      err = DMPlexClearLabelStratum(_dm, "vtk", 2);PYLITH_CHECK_ERROR(err);
    } // if
    err = DMDestroy(&_dm);PYLITH_CHECK_ERROR(err);
  } // if
  _isOpen = false;

  DataWriter<mesh_type, field_type>::close();

  PYLITH_METHOD_END;
} // close

// ----------------------------------------------------------------------
// Prepare file for data at a new time step.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type,field_type>::openTimeStep(const PylithScalar t,
								  const mesh_type& mesh,
								  const char* label,
								  const int labelId)
{ // openTimeStep
  PYLITH_METHOD_BEGIN;

  assert(_dm && _dm == mesh.dmMesh());
  assert(_isOpen && !_isOpenTimeStep);

  PetscErrorCode err = 0;
    
  const std::string& filename = _vtkFilename(t);

  err = PetscViewerCreate(mesh.comm(), &_viewer);PYLITH_CHECK_ERROR(err);
  err = PetscViewerSetType(_viewer, PETSCVIEWERVTK);PYLITH_CHECK_ERROR(err);
  err = PetscViewerSetFormat(_viewer, PETSC_VIEWER_ASCII_VTK);PYLITH_CHECK_ERROR(err);
  err = PetscViewerFileSetName(_viewer, filename.c_str());PYLITH_CHECK_ERROR(err);
  
  // Increment reference count on mesh DM, because the viewer destroys the DM.
  assert(_dm);
  err = PetscObjectReference((PetscObject) _dm);PYLITH_CHECK_ERROR(err);
  
  _isOpenTimeStep = true;

  PYLITH_METHOD_END;
} // openTimeStep

// ----------------------------------------------------------------------
/// Cleanup after writing data for a time step.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type,field_type>::closeTimeStep(void)
{ // closeTimeStep
  PYLITH_METHOD_BEGIN;

  // Account for possibility that no fields were written, so viewer doesn't have handle to DM.
  if (_isOpenTimeStep && !_wroteVertexHeader && !_wroteCellHeader) {
    // No fields written, so must manually dereference the mesh DM.
    PetscErrorCode err = PetscObjectDereference((PetscObject) _dm);PYLITH_CHECK_ERROR(err);
  } // if
  
  PetscErrorCode err = PetscViewerDestroy(&_viewer);PYLITH_CHECK_ERROR(err);
  _isOpenTimeStep = false;
  _wroteVertexHeader = false;
  _wroteCellHeader = false;
  
  PYLITH_METHOD_END;
} // closeTimeStep

// ----------------------------------------------------------------------
// Write field over vertices to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type,field_type>::writeVertexField(const PylithScalar t,
								      field_type& field,
								      const mesh_type& mesh)
{ // writeVertexField
  PYLITH_METHOD_BEGIN;

  assert(_dm && _dm == mesh.dmMesh());
  assert(_isOpen && _isOpenTimeStep);

  // Could check the field.petscSection() matches the default section from VecGetDM().
  Vec v = field.localVector();assert(v);

  // :KLUDGE: MATT You have a note that this is not fully implemented!
  //
  // Will change to just VecView() once I setup the vectors correctly
  // (use VecSetOperation() to change the view method).
  PetscViewerVTKFieldType ft = field.vectorFieldType() != topology::FieldBase::VECTOR ? PETSC_VTK_POINT_FIELD : PETSC_VTK_POINT_VECTOR_FIELD;
  PetscErrorCode err = PetscViewerVTKAddField(_viewer, (PetscObject) _dm, DMPlexVTKWriteAll, ft, (PetscObject) v);PYLITH_CHECK_ERROR(err);
  err = PetscObjectReference((PetscObject) v);PYLITH_CHECK_ERROR(err); /* Needed because viewer destroys the Vec */
  
  _wroteVertexHeader = true;

  PYLITH_METHOD_END;
} // writeVertexField

// ----------------------------------------------------------------------
// Write field over cells to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterVTK<mesh_type,field_type>::writeCellField(const PylithScalar t,
								    field_type& field,
								    const char* label,
								    const int labelId)
{ // writeCellField
  PYLITH_METHOD_BEGIN;

  assert(_dm && _dm == field.mesh().dmMesh());
  assert(_isOpen && _isOpenTimeStep);
  PetscVec v = field.localVector();assert(v);

  // :KLUDGE: MATT You have a note that this is not fully implemented!
  //
  // Will change to just VecView() once I setup the vectors correctly
  // (use VecSetOperation() to change the view).
  PetscViewerVTKFieldType ft = field.vectorFieldType() != topology::FieldBase::VECTOR ? PETSC_VTK_CELL_FIELD : PETSC_VTK_CELL_VECTOR_FIELD;
  PetscErrorCode err = PetscViewerVTKAddField(_viewer, (PetscObject) _dm, DMPlexVTKWriteAll, ft, (PetscObject) v); PYLITH_CHECK_ERROR(err);
  err = PetscObjectReference((PetscObject) v);PYLITH_CHECK_ERROR(err); /* Needed because viewer destroys the Vec */
  
  _wroteCellHeader = true;

  PYLITH_METHOD_END;
} // writeCellField

// ----------------------------------------------------------------------
// Generate filename for VTK file.
template<typename mesh_type, typename field_type>
std::string
pylith::meshio::DataWriterVTK<mesh_type,field_type>::_vtkFilename(const PylithScalar t) const
{ // _vtkFilename
  PYLITH_METHOD_BEGIN;

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

  PYLITH_METHOD_RETURN(std::string(filename.str()));
} // _vtkFilename


// End of file 
