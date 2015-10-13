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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "DataWriterVTK.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // HOLDSA Fields
#include "pylith/topology/Stratum.hh" // USES StratumIS

#include <petscdmplex.h>

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

extern
PetscErrorCode DMPlexVTKWriteAll(PetscObject odm, 
				 PetscViewer viewer);

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::DataWriterVTK::DataWriterVTK(void) :
  _timeConstant(1.0),
  _filename("output.vtk"),
  _timeFormat("%f"),
  _viewer(NULL),
  _dm(NULL),
  _vertexFieldCache(0),
  _cellFieldCache(0),
  _precision(6),
  _isOpen(false),
  _isOpenTimeStep(false),
  _wroteVertexHeader(false),
  _wroteCellHeader(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriterVTK::~DataWriterVTK(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::DataWriterVTK::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  delete _vertexFieldCache; _vertexFieldCache = 0;
  delete _cellFieldCache; _cellFieldCache = 0;

  closeTimeStep(); // Insure time step is closed.
  close(); // Insure clean up.
  DataWriter::deallocate();

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::DataWriterVTK::DataWriterVTK(const DataWriterVTK& w) :
  DataWriter(w),
  _timeConstant(w._timeConstant),
  _filename(w._filename),
  _timeFormat(w._timeFormat),
  _viewer(NULL),
  _dm(NULL),
  _vertexFieldCache(0),
  _cellFieldCache(0),
  _precision(w._precision),
  _isOpen(w._isOpen),
  _isOpenTimeStep(w._isOpenTimeStep),
  _wroteVertexHeader(w._wroteVertexHeader),
  _wroteCellHeader(w._wroteCellHeader)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Set value used to normalize time stamp in name of VTK file.
void
pylith::meshio::DataWriterVTK::timeConstant(const PylithScalar value)
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
void
pylith::meshio::DataWriterVTK::precision(const int value)
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
void
pylith::meshio::DataWriterVTK::open(const topology::Mesh& mesh,
				    const int numTimeSteps,
				    const char* label,
				    const int labelId)
{ // open
  PYLITH_METHOD_BEGIN;

  DataWriter::open(mesh, numTimeSteps, label, labelId);

  // Save handle for actions required in closeTimeStep() and close();
  PetscErrorCode err = 0;
  err = DMDestroy(&_dm);PYLITH_CHECK_ERROR(err);
  _dm = mesh.dmMesh();assert(_dm);
  err = PetscObjectReference((PetscObject) _dm);PYLITH_CHECK_ERROR(err);

  _isOpen = true;

  PYLITH_METHOD_END;
} // open

// ----------------------------------------------------------------------
// Close output files.
void
pylith::meshio::DataWriterVTK::close(void)
{ // close
  PYLITH_METHOD_BEGIN;

  if (_isOpen) {
    assert(_dm);
    PetscErrorCode err = DMDestroy(&_dm);PYLITH_CHECK_ERROR(err);
  } // if

  // Clear field caches.
  if (_vertexFieldCache) {
    _vertexFieldCache->deallocate();
  } // if
  if (_cellFieldCache) {
    _cellFieldCache->deallocate();
  } // if

  DataWriter::close();

  _isOpen = false;

  PYLITH_METHOD_END;
} // close

// ----------------------------------------------------------------------
// Prepare file for data at a new time step.
void
pylith::meshio::DataWriterVTK::openTimeStep(const PylithScalar t,
					    const topology::Mesh& mesh,
					    const char* label,
					    const int labelId)
{ // openTimeStep
  PYLITH_METHOD_BEGIN;

  assert(_dm && _dm == mesh.dmMesh());
  assert(_isOpen && !_isOpenTimeStep);

  PetscErrorCode err = 0;
    
  const std::string& filename = _vtkFilename(t);

  // Create VTK label in DM: Cleared in closeTimeStep().
  if (label) {
    const bool includeOnlyCells = true;
    topology::StratumIS cellsIS(_dm, label, labelId, includeOnlyCells);
    const PetscInt ncells = cellsIS.size();
    const PetscInt* cells = cellsIS.points();
    DMLabel label;

    err = DMPlexCreateLabel(_dm, "vtk");PYLITH_CHECK_ERROR(err);
    err = DMPlexGetLabel(_dm, "vtk", &label);PYLITH_CHECK_ERROR(err);
    err = DMLabelClearStratum(label, 1);PYLITH_CHECK_ERROR(err);
    for (PetscInt c=0; c < ncells; ++c) {
      err = DMLabelSetValue(label, cells[c], 1);PYLITH_CHECK_ERROR(err);
    } // for

  } // if

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
void
pylith::meshio::DataWriterVTK::closeTimeStep(void)
{ // closeTimeStep
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err = 0;

  // Destroy the viewer (which also writes the file).
  err = PetscViewerDestroy(&_viewer);PYLITH_CHECK_ERROR(err);

  // Account for possibility that no fields were written, so viewer doesn't have handle to DM.
  if (_isOpenTimeStep && !_wroteVertexHeader && !_wroteCellHeader) {
    // No fields written, so must manually dereference the mesh DM.
    err = PetscObjectDereference((PetscObject) _dm);PYLITH_CHECK_ERROR(err);
  } // if
  
  // Remove label
  if (_isOpenTimeStep) {
    assert(_dm);
    PetscBool hasLabel = PETSC_FALSE;
    err = DMPlexHasLabel(_dm, "vtk", &hasLabel);PYLITH_CHECK_ERROR(err);
    if (hasLabel) {
      err = DMPlexClearLabelStratum(_dm, "vtk", 1);PYLITH_CHECK_ERROR(err);
      err = DMPlexClearLabelStratum(_dm, "vtk", 2);PYLITH_CHECK_ERROR(err);
    } // if
  } // if

  _isOpenTimeStep = false;
  _wroteVertexHeader = false;
  _wroteCellHeader = false;
  
  PYLITH_METHOD_END;
} // closeTimeStep

// ----------------------------------------------------------------------
// Write field over vertices to file.
void
pylith::meshio::DataWriterVTK::writeVertexField(const PylithScalar t,
						topology::Field& field,
						const topology::Mesh& mesh)
{ // writeVertexField
  PYLITH_METHOD_BEGIN;

  assert(_dm && _dm == mesh.dmMesh());
  assert(_isOpen && _isOpenTimeStep);

  // Cache vertex field since PETSc writer collects all fields and
  // then writes file. Caching the field locally allows the output
  // manager to reuse fields as buffers.
  if (!_vertexFieldCache) {
    _vertexFieldCache = new topology::Fields(field.mesh());assert(_vertexFieldCache);
  } // if/else
  const char* fieldLabel = field.label();
  if (!_vertexFieldCache->hasField(fieldLabel)) {
    _vertexFieldCache->add(fieldLabel, fieldLabel);
    topology::Field& fieldCached = _vertexFieldCache->get(fieldLabel);
    fieldCached.cloneSection(field);
  } // if
  topology::Field& fieldCached = _vertexFieldCache->get(fieldLabel);
  assert(fieldCached.sectionSize() == field.sectionSize());
  fieldCached.copy(field);

  // Could check the field.localSection() matches the default section from VecGetDM().
  PetscVec fieldVec = fieldCached.localVector();assert(fieldVec);

  // :KLUDGE: MATT You have a note that this is not fully implemented!
  //
  // Will change to just VecView() once I setup the vectors correctly
  // (use VecSetOperation() to change the view method).
  PetscViewerVTKFieldType ft = fieldCached.vectorFieldType() != topology::FieldBase::VECTOR ? PETSC_VTK_POINT_FIELD : PETSC_VTK_POINT_VECTOR_FIELD;
  PetscErrorCode err = PetscViewerVTKAddField(_viewer, (PetscObject) _dm, DMPlexVTKWriteAll, ft, (PetscObject) fieldVec);PYLITH_CHECK_ERROR(err);
  err = PetscObjectReference((PetscObject) fieldVec);PYLITH_CHECK_ERROR(err); // Viewer destroys Vec
  
  _wroteVertexHeader = true;

  PYLITH_METHOD_END;
} // writeVertexField

// ----------------------------------------------------------------------
// Write field over cells to file.
void
pylith::meshio::DataWriterVTK::writeCellField(const PylithScalar t,
					      topology::Field& field,
					      const char* label,
					      const int labelId)
{ // writeCellField
  PYLITH_METHOD_BEGIN;

  assert(_dm && _dm == field.mesh().dmMesh());
  assert(_isOpen && _isOpenTimeStep);

  // Cache cell field since PETSc writer collects all fields and
  // then writes file. Caching the field locally allows the output
  // manager to reuse fields as buffers.
  if (!_cellFieldCache) {
    _cellFieldCache = new topology::Fields(field.mesh());assert(_cellFieldCache);
  } // if/else
  const char* fieldLabel = field.label();
  if (!_cellFieldCache->hasField(fieldLabel)) {
    _cellFieldCache->add(fieldLabel, fieldLabel);
    topology::Field& fieldCached = _cellFieldCache->get(fieldLabel);
    fieldCached.cloneSection(field);
  } // if
  topology::Field& fieldCached = _cellFieldCache->get(fieldLabel);
  assert(fieldCached.sectionSize() == field.sectionSize());
  fieldCached.copy(field);

  // Could check the field.localSection() matches the default section from VecGetDM().
  PetscVec fieldVec = fieldCached.localVector();assert(fieldVec);

  // :KLUDGE: MATT You have a note that this is not fully implemented!
  //
  // Will change to just VecView() once I setup the vectors correctly
  // (use VecSetOperation() to change the view).
  PetscViewerVTKFieldType ft = fieldCached.vectorFieldType() != topology::FieldBase::VECTOR ? PETSC_VTK_CELL_FIELD : PETSC_VTK_CELL_VECTOR_FIELD;
  PetscErrorCode err = PetscViewerVTKAddField(_viewer, (PetscObject) _dm, DMPlexVTKWriteAll, ft, (PetscObject) fieldVec); PYLITH_CHECK_ERROR(err);
  err = PetscObjectReference((PetscObject) fieldVec);PYLITH_CHECK_ERROR(err); // Viewer destroys Vec
  
  _wroteCellHeader = true;

  PYLITH_METHOD_END;
} // writeCellField

// ----------------------------------------------------------------------
// Generate filename for VTK file.
std::string
pylith::meshio::DataWriterVTK::_vtkFilename(const PylithScalar t) const
{ // _vtkFilename
  PYLITH_METHOD_BEGIN;

  std::ostringstream filename;
  const int indexExt = _filename.find(".vtk");
  const int numTimeSteps = DataWriter::_numTimeSteps;
  if (numTimeSteps > 0) {
    // If data with multiple time steps, then add time stamp to filename
    char sbuffer[256];
    sprintf(sbuffer, _timeFormat.c_str(), t/_timeConstant);
    std::string timestamp(sbuffer);
    const size_t pos = timestamp.find(".");
    if (pos != std::string::npos) {
      timestamp.erase(pos, 1);
    } // if
    filename << std::string(_filename, 0, indexExt) << "_t" << timestamp << ".vtk";
  } else
    filename << std::string(_filename, 0, indexExt) << "_info.vtk";

  PYLITH_METHOD_RETURN(std::string(filename.str()));
} // _vtkFilename


// End of file 
