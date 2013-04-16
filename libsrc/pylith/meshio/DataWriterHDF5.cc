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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "HDF5.hh" // USES HDF5
#include "Xdmf.hh" // USES Xdmf
#include "petscviewerhdf5.h"

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::DataWriterHDF5(void) :
  _filename("output.h5"),
  _viewer(0),
  _tstamp(0),
  _tstampIndex(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::~DataWriterHDF5(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type, field_type>::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  DataWriter<mesh_type, field_type>::deallocate();

  PetscErrorCode err = 0;
  err = PetscViewerDestroy(&_viewer); CHECK_PETSC_ERROR(err);
  err = VecDestroy(&_tstamp); CHECK_PETSC_ERROR(err);

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor.
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::DataWriterHDF5(const DataWriterHDF5<mesh_type, field_type>& w) :
  DataWriter<mesh_type, field_type>(w),
  _filename(w._filename),
  _viewer(0),
  _tstamp(0),
  _tstampIndex(0)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Prepare file for data at a new time step.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::open(const mesh_type& mesh,
							   const int numTimeSteps,
							   const char* label,
							   const int labelId)
{ // open
  PYLITH_METHOD_BEGIN;

  DataWriter<mesh_type, field_type>::open(mesh, numTimeSteps, label, labelId);

  try {
    PetscErrorCode err = 0;

    deallocate();
    
    const std::string& filename = _hdf5Filename();

    _timesteps.clear();
    _tstampIndex = 0;
    PetscMPIInt commRank;
    err = MPI_Comm_rank(mesh.comm(), &commRank);CHECK_PETSC_ERROR(err);
    const int localSize = (!commRank) ? 1 : 0;
    err = VecCreateMPI(mesh.comm(), localSize, 1, &_tstamp);CHECK_PETSC_ERROR(err);assert(_tstamp);
    err = VecSetBlockSize(_tstamp, 1); CHECK_PETSC_ERROR(err);CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) _tstamp, "time");CHECK_PETSC_ERROR(err);

    err = PetscViewerHDF5Open(mesh.comm(), filename.c_str(), FILE_MODE_WRITE, &_viewer);CHECK_PETSC_ERROR(err);

    const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();assert(cs);

    const char *context = DataWriter<mesh_type, field_type>::_context.c_str();
    PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
    PetscDM dmCoord = NULL;
    PetscVec coordinates = NULL; 
    PetscReal lengthScale;
    topology::FieldBase::Metadata metadata;

    metadata.label = "vertices";
    metadata.vectorFieldType = topology::FieldBase::VECTOR;
    err = DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale);CHECK_PETSC_ERROR(err);
    err = DMGetCoordinateDM(dmMesh, &dmCoord);CHECK_PETSC_ERROR(err);assert(dmCoord);
    err = PetscObjectReference((PetscObject) dmCoord);CHECK_PETSC_ERROR(err);
    err = DMGetCoordinatesLocal(dmMesh, &coordinates);CHECK_PETSC_ERROR(err);
    topology::Field<mesh_type> field(mesh, dmCoord, coordinates, metadata);
    field.createScatterWithBC(mesh, "", 0, metadata.label.c_str());
    field.scatterSectionToVector(metadata.label.c_str());
    PetscVec coordVector = field.vector(metadata.label.c_str());assert(coordVector);
    err = VecScale(coordVector, lengthScale);CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PushGroup(_viewer, "/geometry");CHECK_PETSC_ERROR(err);
    err = VecView(coordVector, _viewer);CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer); CHECK_PETSC_ERROR(err);

    PetscInt vStart, vEnd, cStart, cEnd, cMax, dof, conesSize, numCorners, numCornersLocal = 0;

    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
    err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
    err = DMPlexGetHybridBounds(dmMesh, &cMax, PETSC_NULL, PETSC_NULL, PETSC_NULL);CHECK_PETSC_ERROR(err);
    if (cMax >= 0) {
      cEnd = PetscMin(cEnd, cMax);
    } // if
    for(PetscInt cell = cStart; cell < cEnd; ++cell) {
      PetscInt *closure = NULL;
      PetscInt  closureSize, v;

      err = DMPlexGetTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
      numCornersLocal = 0;
      for (v = 0; v < closureSize*2; v += 2) {
        if ((closure[v] >= vStart) && (closure[v] < vEnd)) {
          ++numCornersLocal;
        } // if
      } // for
      err = DMPlexRestoreTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
      if (numCornersLocal)
	break;
    } // for
    err = MPI_Allreduce(&numCornersLocal, &numCorners, 1, MPIU_INT, MPI_MAX, mesh.comm());CHECK_PETSC_ERROR(err);

    if (label) {
      conesSize = 0;
      for(PetscInt cell = cStart; cell < cEnd; ++cell) {
        PetscInt value;
        err = DMPlexGetLabelValue(dmMesh, label, cell, &value);CHECK_PETSC_ERROR(err);
        if (value == labelId)
	  ++conesSize;
      } // for
      conesSize *= numCorners;
    } else {
      conesSize = (cEnd - cStart)*numCorners;
    } // if/else
    CHKMEMA; // MATT Why is this here?

    PetscIS globalVertexNumbers = NULL;
    const PetscInt *gvertex = NULL;
    PetscVec cellVec = NULL;
    PetscScalar *vertices = NULL;

    err = DMPlexGetVertexNumbering(dmMesh, &globalVertexNumbers);CHECK_PETSC_ERROR(err);
    err = ISGetIndices(globalVertexNumbers, &gvertex);CHECK_PETSC_ERROR(err);
    err = VecCreate(mesh.comm(), &cellVec);CHECK_PETSC_ERROR(err);
    err = VecSetSizes(cellVec, conesSize, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(cellVec, numCorners);CHECK_PETSC_ERROR(err);
    err = VecSetFromOptions(cellVec);CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) cellVec, "cells");CHECK_PETSC_ERROR(err);
    err = VecGetArray(cellVec, &vertices);CHECK_PETSC_ERROR(err);
    for(PetscInt cell = cStart, v = 0; cell < cEnd; ++cell) {
      PetscInt *closure = NULL;
      PetscInt  closureSize, p;

      if (label) {
        PetscInt value;
        err = DMPlexGetLabelValue(dmMesh, label, cell, &value);CHECK_PETSC_ERROR(err);
        if (value != labelId)
	  continue;
      } // if

      err = DMPlexGetTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
      for(p = 0; p < closureSize*2; p += 2) {
        if ((closure[p] >= vStart) && (closure[p] < vEnd)) {
          const PetscInt gv = gvertex[closure[p] - vStart];
          vertices[v++] = gv < 0 ? -(gv+1) : gv;
        } // if
      } // for
      err = DMPlexRestoreTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
      //assert(v == (cell-cStart+1)*numCorners); :TODO: MATT Why is this here?
    } // for
    CHKMEMA; // MATT why is this here?
    err = VecRestoreArray(cellVec, &vertices);CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PushGroup(_viewer, "/topology");CHECK_PETSC_ERROR(err);
    err = VecView(cellVec, _viewer);CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer);CHECK_PETSC_ERROR(err);
    err = VecDestroy(&cellVec);CHECK_PETSC_ERROR(err);
    err = ISRestoreIndices(globalVertexNumbers, &gvertex);CHECK_PETSC_ERROR(err);

    hid_t h5 = -1;
    err = PetscViewerHDF5GetFileId(_viewer, &h5); CHECK_PETSC_ERROR(err);
    assert(h5 >= 0);
    const int cellDim = mesh.dimension();
    HDF5::writeAttribute(h5, "/topology/cells", "cell_dim", (void*)&cellDim, H5T_NATIVE_INT);
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while opening HDF5 file " << _hdf5Filename() << ".\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Unknown error while opening HDF5 file " << _hdf5Filename() << ".\n";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_END;
} // open

// ----------------------------------------------------------------------
// Close output files.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::close(void)
{ // close
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err = 0;
  err = PetscViewerDestroy(&_viewer); CHECK_PETSC_ERROR(err);
  err = VecDestroy(&_tstamp); CHECK_PETSC_ERROR(err);assert(!_tstamp);

  _timesteps.clear();
  _tstampIndex = 0;

  int rank = 0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    Xdmf metafile;
    const std::string& hdf5filename = _hdf5Filename();
    const int indexExt = hdf5filename.find(".h5");
    std::string xdmfFilename = std::string(hdf5filename, 0, indexExt) + ".xmf";
    metafile.write(xdmfFilename.c_str(), _hdf5Filename().c_str());
  } // if

  PYLITH_METHOD_END;
} // close

// ----------------------------------------------------------------------
// Write field over vertices to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::writeVertexField(const PylithScalar t,
								       field_type& field,
								       const mesh_type& mesh)
{ // writeVertexField
  PYLITH_METHOD_BEGIN;

  assert(_viewer);

  try {
    PetscErrorCode err;

    const char* context  = DataWriter<mesh_type, field_type>::_context.c_str();

    field.createScatterWithBC(mesh, "", 0, context);
    field.scatterSectionToVector(context);
    PetscVec vector = field.vector(context);assert(vector);

    if (_timesteps.find(field.label()) == _timesteps.end())
      _timesteps[field.label()] = 0;
    else
      _timesteps[field.label()] += 1;
    const int istep = _timesteps[field.label()];
    // Add time stamp to "/time" if necessary.
    PetscMPIInt commRank;
    err = MPI_Comm_rank(mesh.comm(), &commRank);CHECK_PETSC_ERROR(err);
    if (_tstampIndex == istep)
      _writeTimeStamp(t, commRank);

#if 0 // :TODO: MATT What is this doing here?
    const int spaceDim = mesh.coordsys()->spaceDim();
    PetscInt  bs;
    err = VecGetBlockSize(vector, &bs);CHECK_PETSC_ERROR(err);
    switch (field.vectorFieldType()) {
    case pylith::topology::FieldBase::VECTOR:
      if (bs % spaceDim) {
	CHECK_PETSC_ERROR(PETSC_ERR_ARG_WRONG);
      } // if
      break;
    case pylith::topology::FieldBase::TENSOR:
      if (bs % spaceDim) {
      CHECK_PETSC_ERROR(PETSC_ERR_ARG_WRONG);
      } // if
      break;
    default:
      if (bs > 1) {
	CHECK_PETSC_ERROR(PETSC_ERR_ARG_WRONG);
      } // if
      break;
    } // switch
#endif

    err = PetscViewerHDF5PushGroup(_viewer, "/vertex_fields");CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5SetTimestep(_viewer, istep);CHECK_PETSC_ERROR(err);
    err = VecView(vector, _viewer);CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer);CHECK_PETSC_ERROR(err);

    if (0 == istep) {
      hid_t h5 = -1;
      err = PetscViewerHDF5GetFileId(_viewer, &h5); CHECK_PETSC_ERROR(err);
      assert(h5 >= 0);
      std::string fullName = std::string("/vertex_fields/") + field.label();
      HDF5::writeAttribute(h5, fullName.c_str(), "vector_field_type",
			   topology::FieldBase::vectorFieldString(field.vectorFieldType()));
    } // if

  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to HDF5 file '" << _hdf5Filename() << "'.\n" << err.what();
    throw std::runtime_error(msg.str());

  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to HDF5 file '" << _hdf5Filename() << "'.\n";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_END;
} // writeVertexField

// ----------------------------------------------------------------------
// Write field over cells to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::writeCellField(const PylithScalar t,
								     field_type& field,
								     const char* label,
								     const int labelId)
{ // writeCellField
  PYLITH_METHOD_BEGIN;

  assert(_viewer);
  
  try {
    const char* context = DataWriter<mesh_type, field_type>::_context.c_str();
    PetscErrorCode err = 0;

    field.createScatterWithBC(field.mesh(), label ? label : "", labelId, context);
    field.scatterSectionToVector(context);
    PetscVec vector = field.vector(context);assert(vector);

    if (_timesteps.find(field.label()) == _timesteps.end())
      _timesteps[field.label()] = 0;
    else
      _timesteps[field.label()] += 1;
    const int istep = _timesteps[field.label()];
    // Add time stamp to "/time" if necessary.
    PetscMPIInt commRank;
    err = MPI_Comm_rank(field.mesh().comm(), &commRank);CHECK_PETSC_ERROR(err);
    if (_tstampIndex == istep)
      _writeTimeStamp(t, commRank);

    err = PetscViewerHDF5PushGroup(_viewer, "/cell_fields");CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5SetTimestep(_viewer, istep);CHECK_PETSC_ERROR(err);
    err = VecView(vector, _viewer);CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer);CHECK_PETSC_ERROR(err);

    if (0 == istep) {
      hid_t h5 = -1;
      err = PetscViewerHDF5GetFileId(_viewer, &h5); CHECK_PETSC_ERROR(err);
      assert(h5 >= 0);
      std::string fullName = std::string("/cell_fields/") + field.label();
      HDF5::writeAttribute(h5, fullName.c_str(), "vector_field_type",
			   topology::FieldBase::vectorFieldString(field.vectorFieldType()));
    } // if
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to HDF5 file '" << _hdf5Filename() << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to HDF5 file '" << _hdf5Filename() << "'.\n";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_END;
} // writeCellField

// ----------------------------------------------------------------------
// Generate filename for HDF5 file.
template<typename mesh_type, typename field_type>
std::string
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::_hdf5Filename(void) const
{ // _hdf5Filename
  PYLITH_METHOD_BEGIN;

  std::ostringstream filename;
  const int indexExt = _filename.find(".h5");
  const int numTimeSteps = DataWriter<mesh_type, field_type>::_numTimeSteps;
  if (0 == numTimeSteps) {
    filename << std::string(_filename, 0, indexExt) << "_info.h5";
  } else {
    filename << _filename;
  } // if/else

  PYLITH_METHOD_RETURN(std::string(filename.str()));
} // _hdf5Filename


// ----------------------------------------------------------------------
// Write time stamp to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::_writeTimeStamp(const PylithScalar t,
								      const int commRank)
{ // _writeTimeStamp
  assert(_tstamp);
  PetscErrorCode err = 0;

  if (0 == commRank) {
    const PylithScalar tDim = t * DataWriter<mesh_type, field_type>::_timeScale;
    err = VecSetValue(_tstamp, 0, tDim, INSERT_VALUES); CHECK_PETSC_ERROR(err);
  } // if
  err = VecAssemblyBegin(_tstamp); CHECK_PETSC_ERROR(err);
  err = VecAssemblyEnd(_tstamp); CHECK_PETSC_ERROR(err);
  
  err = PetscViewerHDF5PushGroup(_viewer, "/"); CHECK_PETSC_ERROR(err);
  err = PetscViewerHDF5SetTimestep(_viewer, _tstampIndex); CHECK_PETSC_ERROR(err);
  err = VecView(_tstamp, _viewer); CHECK_PETSC_ERROR(err);
  err = PetscViewerHDF5PopGroup(_viewer); CHECK_PETSC_ERROR(err);

  _tstampIndex++;
} // _writeTimeStamp


// End of file 
