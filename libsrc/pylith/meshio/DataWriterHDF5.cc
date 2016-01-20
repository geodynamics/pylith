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

#include "DataWriterHDF5.hh" // Implementation of class methods

#include "HDF5.hh" // USES HDF5
#include "Xdmf.hh" // USES Xdmf

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "petscviewerhdf5.h"

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

extern "C" {
  extern PetscErrorCode VecView_Seq(Vec, PetscViewer);
  extern PetscErrorCode VecView_MPI(Vec, PetscViewer);
}

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::DataWriterHDF5::DataWriterHDF5(void) :
  _filename("output.h5"),
  _viewer(0),
  _tstamp(0),
  _tstampIndex(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriterHDF5::~DataWriterHDF5(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::DataWriterHDF5::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  DataWriter::deallocate();

  PetscErrorCode err = 0;
  err = PetscViewerDestroy(&_viewer); PYLITH_CHECK_ERROR(err);assert(!_viewer);
  err = VecDestroy(&_tstamp); PYLITH_CHECK_ERROR(err);assert(!_tstamp);

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::DataWriterHDF5::DataWriterHDF5(const DataWriterHDF5& w) :
  DataWriter(w),
  _filename(w._filename),
  _viewer(0),
  _tstamp(0),
  _tstampIndex(0)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Prepare file for data at a new time step.
void
pylith::meshio::DataWriterHDF5::open(const topology::Mesh& mesh,
				     const int numTimeSteps,
				     const char* label,
				     const int labelId)
{ // open
  PYLITH_METHOD_BEGIN;

  DataWriter::open(mesh, numTimeSteps, label, labelId);

  try {
    PetscErrorCode err = 0;

    deallocate();
    
    const std::string& filename = _hdf5Filename();

    _timesteps.clear();
    _tstampIndex = 0;
    PetscMPIInt commRank;
    err = MPI_Comm_rank(mesh.comm(), &commRank);PYLITH_CHECK_ERROR(err);
    const int localSize = (!commRank) ? 1 : 0;
    err = VecCreateMPI(mesh.comm(), localSize, 1, &_tstamp);PYLITH_CHECK_ERROR(err);assert(_tstamp);
    err = VecSetBlockSize(_tstamp, 1); PYLITH_CHECK_ERROR(err);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) _tstamp, "time");PYLITH_CHECK_ERROR(err);

    err = PetscViewerHDF5Open(mesh.comm(), filename.c_str(), FILE_MODE_WRITE, &_viewer);PYLITH_CHECK_ERROR(err);
    err = PetscViewerHDF5SetBaseDimension2(_viewer, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
    
    const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();assert(cs);

    PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
    PetscDM dmCoord = NULL;
    PetscVec coordinates = NULL; 
    PetscReal lengthScale;
    topology::FieldBase::Metadata metadata;

    metadata.label = "vertices";
    metadata.vectorFieldType = topology::FieldBase::VECTOR;
    err = DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMGetCoordinateDM(dmMesh, &dmCoord);PYLITH_CHECK_ERROR(err);assert(dmCoord);
    err = PetscObjectReference((PetscObject) dmCoord);PYLITH_CHECK_ERROR(err);
    err = DMGetCoordinatesLocal(dmMesh, &coordinates);PYLITH_CHECK_ERROR(err);
    topology::Field coordinatesField(mesh, dmCoord, coordinates, metadata);
    coordinatesField.createScatterWithBC(mesh, "", 0, metadata.label.c_str());
    coordinatesField.scatterLocalToGlobal(metadata.label.c_str());
    PetscVec coordVector = coordinatesField.vector(metadata.label.c_str());assert(coordVector);
    err = VecScale(coordVector, lengthScale);PYLITH_CHECK_ERROR(err);
    err = PetscViewerHDF5PushGroup(_viewer, "/geometry");PYLITH_CHECK_ERROR(err);
#if 0
    err = VecView(coordVector, _viewer);PYLITH_CHECK_ERROR(err);
#else
    PetscBool isseq;
    err = PetscObjectTypeCompare((PetscObject) coordVector, VECSEQ, &isseq);PYLITH_CHECK_ERROR(err);
    if (isseq) {err = VecView_Seq(coordVector, _viewer);PYLITH_CHECK_ERROR(err);}
    else       {err = VecView_MPI(coordVector, _viewer);PYLITH_CHECK_ERROR(err);}
#endif
    err = PetscViewerHDF5PopGroup(_viewer); PYLITH_CHECK_ERROR(err);

    PetscInt vStart, vEnd, cellHeight, cStart, cEnd, cMax, conesSize, numCorners, numCornersLocal = 0;

    err = DMPlexGetVTKCellHeight(dmMesh, &cellHeight);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetHeightStratum(dmMesh, cellHeight, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetHybridBounds(dmMesh, &cMax, PETSC_NULL, PETSC_NULL, PETSC_NULL);PYLITH_CHECK_ERROR(err);
    if (cMax >= 0) {
      cEnd = PetscMin(cEnd, cMax);
    } // if
    for(PetscInt cell = cStart; cell < cEnd; ++cell) {
      PetscInt *closure = NULL;
      PetscInt  closureSize, v;

      err = DMPlexGetTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
      numCornersLocal = 0;
      for (v = 0; v < closureSize*2; v += 2) {
        if ((closure[v] >= vStart) && (closure[v] < vEnd)) {
          ++numCornersLocal;
        } // if
      } // for
      err = DMPlexRestoreTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
      if (numCornersLocal)
	break;
    } // for
    err = MPI_Allreduce(&numCornersLocal, &numCorners, 1, MPIU_INT, MPI_MAX, mesh.comm());PYLITH_CHECK_ERROR(err);

    if (label) {
      conesSize = 0;
      for(PetscInt cell = cStart; cell < cEnd; ++cell) {
        PetscInt value;
        err = DMGetLabelValue(dmMesh, label, cell, &value);PYLITH_CHECK_ERROR(err);
        if (value == labelId)
	  ++conesSize;
      } // for
      conesSize *= numCorners;
    } else {
      conesSize = (cEnd - cStart)*numCorners;
    } // if/else

    PetscIS globalVertexNumbers = NULL;
    const PetscInt *gvertex = NULL;
    PetscVec cellVec = NULL;
    PetscScalar *vertices = NULL;
    const PetscInt dim = mesh.dimension();

    err = DMPlexGetVertexNumbering(dmMesh, &globalVertexNumbers);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(globalVertexNumbers, &gvertex);PYLITH_CHECK_ERROR(err);
    err = VecCreate(mesh.comm(), &cellVec);PYLITH_CHECK_ERROR(err);
    err = VecSetSizes(cellVec, conesSize, PETSC_DETERMINE);PYLITH_CHECK_ERROR(err);
    err = VecSetBlockSize(cellVec, numCorners);PYLITH_CHECK_ERROR(err);
    err = VecSetFromOptions(cellVec);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) cellVec, "cells");PYLITH_CHECK_ERROR(err);
    err = VecGetArray(cellVec, &vertices);PYLITH_CHECK_ERROR(err);
    for(PetscInt cell = cStart, v = 0; cell < cEnd; ++cell) {
      PetscInt *closure = NULL;
      PetscInt  closureSize, nC = 0, p;

      if (label) {
        PetscInt value;
        err = DMGetLabelValue(dmMesh, label, cell, &value);PYLITH_CHECK_ERROR(err);
        if (value != labelId) continue;
      } // if

      err = DMPlexGetTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
      for(p = 0; p < closureSize*2; p += 2) {
        if ((closure[p] >= vStart) && (closure[p] < vEnd)) {
          closure[nC++] = closure[p];
        } // if
      } // for
      err = DMPlexInvertCell(dim, nC, closure);PYLITH_CHECK_ERROR(err);
      for (p = 0; p < nC; ++p) {
        const PetscInt gv = gvertex[closure[p] - vStart];
        vertices[v++] = gv < 0 ? -(gv+1) : gv;
      }
      err = DMPlexRestoreTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
      //assert(v == (cell-cStart+1)*numCorners); Would be true without the label check
    } // for
    err = VecRestoreArray(cellVec, &vertices);PYLITH_CHECK_ERROR(err);
    err = PetscViewerHDF5PushGroup(_viewer, "/topology");PYLITH_CHECK_ERROR(err);
    err = VecView(cellVec, _viewer);PYLITH_CHECK_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&cellVec);PYLITH_CHECK_ERROR(err);
    err = ISRestoreIndices(globalVertexNumbers, &gvertex);PYLITH_CHECK_ERROR(err);

    hid_t h5 = -1;
    err = PetscViewerHDF5GetFileId(_viewer, &h5); PYLITH_CHECK_ERROR(err);
    assert(h5 >= 0);
    const int cellDim = mesh.dimension();
    HDF5::writeAttribute(h5, "/topology/cells", "cell_dim", (void*)&cellDim, H5T_NATIVE_INT);
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while opening HDF5 file " << _hdf5Filename() << ".\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Unknown error while opening HDF5 file " << _hdf5Filename() << ".";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_END;
} // open

// ----------------------------------------------------------------------
// Close output files.
void
pylith::meshio::DataWriterHDF5::close(void)
{ // close
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err = 0;
  err = PetscViewerDestroy(&_viewer); PYLITH_CHECK_ERROR(err);assert(!_viewer);
  err = VecDestroy(&_tstamp); PYLITH_CHECK_ERROR(err);assert(!_tstamp);

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
void
pylith::meshio::DataWriterHDF5::writeVertexField(const PylithScalar t,
						 topology::Field& field,
						 const topology::Mesh& mesh)
{ // writeVertexField
  PYLITH_METHOD_BEGIN;

  assert(_viewer);

  try {
    PetscErrorCode err;

    const char* context  = DataWriter::_context.c_str();

    field.createScatterWithBC(mesh, "", 0, context);
    field.scatterLocalToGlobal(context);
    PetscVec vector = field.vector(context);assert(vector);

    if (_timesteps.find(field.label()) == _timesteps.end())
      _timesteps[field.label()] = 0;
    else
      _timesteps[field.label()] += 1;
    const int istep = _timesteps[field.label()];
    // Add time stamp to "/time" if necessary.
    PetscMPIInt commRank;
    err = MPI_Comm_rank(mesh.comm(), &commRank);PYLITH_CHECK_ERROR(err);
    if (_tstampIndex == istep)
      _writeTimeStamp(t, commRank);

    err = PetscViewerHDF5PushGroup(_viewer, "/vertex_fields");PYLITH_CHECK_ERROR(err);
    err = PetscViewerHDF5SetTimestep(_viewer, istep);PYLITH_CHECK_ERROR(err);
#if 0
    err = VecView(vector, _viewer);PYLITH_CHECK_ERROR(err);
#else
    PetscBool isseq;
    err = PetscObjectTypeCompare((PetscObject) vector, VECSEQ, &isseq);PYLITH_CHECK_ERROR(err);
    if (isseq) {err = VecView_Seq(vector, _viewer);PYLITH_CHECK_ERROR(err);}
    else       {err = VecView_MPI(vector, _viewer);PYLITH_CHECK_ERROR(err);}
#endif
    err = PetscViewerHDF5PopGroup(_viewer);PYLITH_CHECK_ERROR(err);

    if (0 == istep) {
      hid_t h5 = -1;
      err = PetscViewerHDF5GetFileId(_viewer, &h5); PYLITH_CHECK_ERROR(err);
      assert(h5 >= 0);
      std::string fullName = std::string("/vertex_fields/") + field.label();
      const char* sattr = topology::FieldBase::vectorFieldString(field.vectorFieldType());
      HDF5::writeAttribute(h5, fullName.c_str(), "vector_field_type", sattr);
    } // if

  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to HDF5 file '" << _hdf5Filename() << "'.\n" << err.what();
    throw std::runtime_error(msg.str());

  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to HDF5 file '" << _hdf5Filename() << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_END;
} // writeVertexField

// ----------------------------------------------------------------------
// Write field over cells to file.
void
pylith::meshio::DataWriterHDF5::writeCellField(const PylithScalar t,
					       topology::Field& field,
					       const char* label,
					       const int labelId)
{ // writeCellField
  PYLITH_METHOD_BEGIN;

  assert(_viewer);
  
  try {
    const char* context = DataWriter::_context.c_str();
    PetscErrorCode err = 0;

    field.createScatterWithBC(field.mesh(), label ? label : "", labelId, context);
    field.scatterLocalToGlobal(context);
    PetscVec vector = field.vector(context);assert(vector);

    if (_timesteps.find(field.label()) == _timesteps.end())
      _timesteps[field.label()] = 0;
    else
      _timesteps[field.label()] += 1;
    const int istep = _timesteps[field.label()];
    // Add time stamp to "/time" if necessary.
    PetscMPIInt commRank;
    err = MPI_Comm_rank(field.mesh().comm(), &commRank);PYLITH_CHECK_ERROR(err);
    if (_tstampIndex == istep)
      _writeTimeStamp(t, commRank);

    err = PetscViewerHDF5PushGroup(_viewer, "/cell_fields");PYLITH_CHECK_ERROR(err);
    err = PetscViewerHDF5SetTimestep(_viewer, istep);PYLITH_CHECK_ERROR(err);
#if 0
    err = VecView(vector, _viewer);PYLITH_CHECK_ERROR(err);
#else
    PetscBool isseq;
    err = PetscObjectTypeCompare((PetscObject) vector, VECSEQ, &isseq);PYLITH_CHECK_ERROR(err);
    if (isseq) {err = VecView_Seq(vector, _viewer);PYLITH_CHECK_ERROR(err);}
    else       {err = VecView_MPI(vector, _viewer);PYLITH_CHECK_ERROR(err);}
#endif
    err = PetscViewerHDF5PopGroup(_viewer);PYLITH_CHECK_ERROR(err);

    if (0 == istep) {
      hid_t h5 = -1;
      err = PetscViewerHDF5GetFileId(_viewer, &h5); PYLITH_CHECK_ERROR(err);
      assert(h5 >= 0);
      std::string fullName = std::string("/cell_fields/") + field.label();
      const char* sattr = topology::FieldBase::vectorFieldString(field.vectorFieldType());
      HDF5::writeAttribute(h5, fullName.c_str(), "vector_field_type", sattr);
    } // if
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to HDF5 file '" << _hdf5Filename() << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to HDF5 file '" << _hdf5Filename() << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_END;
} // writeCellField

// ----------------------------------------------------------------------
// Write dataset with names of points to file.
void
pylith::meshio::DataWriterHDF5::writePointNames(const char* const* names,
						const int numNames,
						const topology::Mesh& mesh)
{ // writePointNames
  PYLITH_METHOD_BEGIN;

  assert(_viewer);

  try {
    hid_t h5 = -1;
    PetscErrorCode err = PetscViewerHDF5GetFileId(_viewer, &h5); PYLITH_CHECK_ERROR(err);
    assert(h5 >= 0);
    HDF5::writeDataset(h5, "/", "stations", names, numNames);
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing stations to HDF5 file '" << _hdf5Filename() << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing stations to HDF5 file '" << _hdf5Filename() << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_END;
} // writePointNames

// ----------------------------------------------------------------------
// Generate filename for HDF5 file.
std::string
pylith::meshio::DataWriterHDF5::_hdf5Filename(void) const
{ // _hdf5Filename
  PYLITH_METHOD_BEGIN;

  std::ostringstream filename;
  const int indexExt = _filename.find(".h5");
  const int numTimeSteps = DataWriter::_numTimeSteps;
  if (0 == numTimeSteps) {
    filename << std::string(_filename, 0, indexExt) << "_info.h5";
  } else {
    filename << _filename;
  } // if/else

  PYLITH_METHOD_RETURN(std::string(filename.str()));
} // _hdf5Filename


// ----------------------------------------------------------------------
// Write time stamp to file.
void
pylith::meshio::DataWriterHDF5::_writeTimeStamp(const PylithScalar t,
						const int commRank)
{ // _writeTimeStamp
  assert(_tstamp);
  PetscErrorCode err = 0;

  if (!commRank) {
    const PylithScalar tDim = t * DataWriter::_timeScale;
    err = VecSetValue(_tstamp, 0, tDim, INSERT_VALUES); PYLITH_CHECK_ERROR(err);
  } // if
  err = VecAssemblyBegin(_tstamp); PYLITH_CHECK_ERROR(err);
  err = VecAssemblyEnd(_tstamp); PYLITH_CHECK_ERROR(err);
  
  err = PetscViewerHDF5PushGroup(_viewer, "/"); PYLITH_CHECK_ERROR(err);
  err = PetscViewerHDF5SetTimestep(_viewer, _tstampIndex); PYLITH_CHECK_ERROR(err);
  err = VecView(_tstamp, _viewer); PYLITH_CHECK_ERROR(err);
  err = PetscViewerHDF5PopGroup(_viewer); PYLITH_CHECK_ERROR(err);

  _tstampIndex++;
} // _writeTimeStamp


// End of file 
