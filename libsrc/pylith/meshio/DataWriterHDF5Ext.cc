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

#include "DataWriterHDF5Ext.hh" // Implementation of class methods

#include "HDF5.hh" // USES HDF5
#include "Xdmf.hh" // USES Xdmf

#include "pylith/topology/Mesh.hh" /// USES Mesh
#include "pylith/topology/Field.hh" /// USES Field
#include "pylith/topology/Stratum.hh" /// USES StratumIS

#include "spatialdata/geocoords/CoordSys.hh" /// USES CoordSys

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

extern "C" {
extern PetscErrorCode VecView_Seq(Vec, PetscViewer);
extern PetscErrorCode VecView_MPI(Vec, PetscViewer);
}

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::DataWriterHDF5Ext::DataWriterHDF5Ext(void) :
  _filename("output.h5"),
  _h5(new HDF5),
  _tstampIndex(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriterHDF5Ext::~DataWriterHDF5Ext(void)
{ // destructor
  delete _h5; _h5 = 0;
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::DataWriterHDF5Ext::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  DataWriter::deallocate();

  PetscErrorCode err = 0;
  const dataset_type::const_iterator& dEnd = _datasets.end();
  for (dataset_type::iterator d_iter=_datasets.begin();
       d_iter != dEnd;
       ++d_iter) {
    err = PetscViewerDestroy(&d_iter->second.viewer);PYLITH_CHECK_ERROR(err);
  } // for

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::DataWriterHDF5Ext::DataWriterHDF5Ext(const DataWriterHDF5Ext& w) :
  DataWriter(w),
  _filename(w._filename),
  _h5(new HDF5),
  _tstampIndex(0)
{ // copy constructor
} // copy constructor

#include <iostream>
// ----------------------------------------------------------------------
// Prepare for writing files.
void
pylith::meshio::DataWriterHDF5Ext::open(const topology::Mesh& mesh,
					const int numTimeSteps,
					const char* label,
					const int labelId)
{ // open
  PYLITH_METHOD_BEGIN;

  assert(_h5);
  _datasets.clear();

  try {
    DataWriter::open(mesh, numTimeSteps, label, labelId);

    PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
    MPI_Comm comm;
    PetscMPIInt commRank;
    PetscErrorCode err = PetscObjectGetComm((PetscObject) dmMesh, &comm);PYLITH_CHECK_ERROR(err);

    err = MPI_Comm_rank(comm, &commRank);PYLITH_CHECK_ERROR(err);
    if (!commRank) {
      _h5->open(_hdf5Filename().c_str(), H5F_ACC_TRUNC);

      // Create groups
      _h5->createGroup("/topology");
      _h5->createGroup("/geometry");
    } // if
    _tstampIndex = 0;

    PetscViewer binaryViewer;
    
    const hid_t scalartype = (sizeof(double) == sizeof(PylithScalar)) ? H5T_IEEE_F64BE : H5T_IEEE_F32BE;

    // Write vertex coordinates
    const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();assert(cs);

    /* TODO Get rid of this and use the createScatterWithBC(numbering) code */
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

    const std::string& filenameVertices = _datasetFilename("vertices");
    err = PetscViewerBinaryOpen(comm, filenameVertices.c_str(), FILE_MODE_WRITE, &binaryViewer);PYLITH_CHECK_ERROR(err);
    err = PetscViewerBinarySetSkipHeader(binaryViewer, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
#if 0
    err = VecView(coordVector, binaryViewer); PYLITH_CHECK_ERROR(err);
#else
    PetscBool isseq;
    err = PetscObjectTypeCompare((PetscObject) coordVector, VECSEQ, &isseq);PYLITH_CHECK_ERROR(err);
    if (isseq) {err = VecView_Seq(coordVector, binaryViewer);PYLITH_CHECK_ERROR(err);}
    else       {err = VecView_MPI(coordVector, binaryViewer);PYLITH_CHECK_ERROR(err);}
#endif
    err = PetscViewerDestroy(&binaryViewer); PYLITH_CHECK_ERROR(err);

    PetscInt vStart, vEnd;
    PetscInt n, numVerticesLocal = 0, numVertices;
    PetscIS globalVertexNumbers = NULL;
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetVertexNumbering(dmMesh, &globalVertexNumbers);PYLITH_CHECK_ERROR(err);
    err = ISGetLocalSize(globalVertexNumbers, &n);PYLITH_CHECK_ERROR(err);
    if (n > 0) {
      const PetscInt *indices = NULL;
      err = ISGetIndices(globalVertexNumbers, &indices);PYLITH_CHECK_ERROR(err);
      for(PetscInt v = 0; v < n; ++v) {
	if (indices[v] >= 0) ++numVerticesLocal;
      } // for
      err = ISRestoreIndices(globalVertexNumbers, &indices);PYLITH_CHECK_ERROR(err);
    } // if
    err = MPI_Allreduce(&numVerticesLocal, &numVertices, 1, MPI_INT, MPI_SUM, comm);PYLITH_CHECK_ERROR(err);
    assert(numVertices > 0);
    
    // Create external dataset for coordinates    
    if (!commRank) {
      const hsize_t ndims = 2;
      hsize_t dims[ndims];
      dims[0] = numVertices;
      dims[1] = cs->spaceDim();
      _h5->createDatasetRawExternal("/geometry", "vertices", filenameVertices.c_str(), dims, ndims, scalartype);
    } // if
    
    // Write cells

    // Account for censored cells
    PetscInt cellHeight, cStart, cEnd, cMax, conesSize, numCells, numCorners, numCornersLocal = 0;
    err = DMPlexGetVTKCellHeight(dmMesh, &cellHeight);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetHeightStratum(dmMesh, cellHeight, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetHybridBounds(dmMesh, &cMax, NULL, NULL, NULL);PYLITH_CHECK_ERROR(err);
    if (cMax >= 0) {
      cEnd = PetscMin(cEnd, cMax);
    } // if
    for(PetscInt cell = cStart; cell < cEnd; ++cell) {
      PetscInt *closure = NULL;
      PetscInt closureSize, v;

      err = DMPlexGetTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
      numCornersLocal = 0;
      for (v = 0; v < closureSize*2; v += 2) {
        if ((closure[v] >= vStart) && (closure[v] < vEnd)) {
          ++numCornersLocal;
        } // if
      } // for
      err = DMPlexRestoreTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
      if (numCornersLocal) break;
    } // for
    err = MPI_Allreduce(&numCornersLocal, &numCorners, 1, MPIU_INT, MPI_MAX, comm);PYLITH_CHECK_ERROR(err);
    if (label) {
      conesSize = 0;
      for(PetscInt cell = cStart; cell < cEnd; ++cell) {
        PetscInt value;

        err = DMPlexGetLabelValue(dmMesh, label, cell, &value);PYLITH_CHECK_ERROR(err);
        if (value == labelId) ++conesSize;
      } // for
      conesSize *= numCorners;
    } else {
      conesSize = (cEnd - cStart)*numCorners;
    } // if/else

    globalVertexNumbers = NULL;
    const PetscInt *gvertex = NULL;
    PetscVec cellVec = NULL;
    PetscScalar *vertices = NULL;
    const PetscInt meshDim = mesh.dimension();

    err = DMPlexGetVertexNumbering(dmMesh, &globalVertexNumbers);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(globalVertexNumbers, &gvertex);PYLITH_CHECK_ERROR(err);
    err = VecCreate(comm, &cellVec);PYLITH_CHECK_ERROR(err);
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

        err = DMPlexGetLabelValue(dmMesh, label, cell, &value);PYLITH_CHECK_ERROR(err);
        if (value != labelId) continue;
      } // if
      err = DMPlexGetTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
      for(p = 0; p < closureSize*2; p += 2) {
        if ((closure[p] >= vStart) && (closure[p] < vEnd)) {
          closure[nC++] = closure[p];
        } // if
      } // for
      err = DMPlexInvertCell(meshDim, nC, closure);PYLITH_CHECK_ERROR(err);
      for (p = 0; p < nC; ++p) {
        const PetscInt gv = gvertex[closure[p] - vStart];
        vertices[v++] = gv < 0 ? -(gv+1) : gv;
      }
      err = DMPlexRestoreTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
      //assert(v == (cell-cStart+1)*numCorners); Would be true without label check
    } // for
    err = ISRestoreIndices(globalVertexNumbers, &gvertex);PYLITH_CHECK_ERROR(err);
    err = VecRestoreArray(cellVec, &vertices);PYLITH_CHECK_ERROR(err);
    err = VecGetSize(cellVec, &numCells);PYLITH_CHECK_ERROR(err);
    numCells /= numCorners;

    const std::string& filenameCells = _datasetFilename("cells");
    err = PetscViewerBinaryOpen(comm, filenameCells.c_str(), FILE_MODE_WRITE, &binaryViewer);PYLITH_CHECK_ERROR(err);
    err = PetscViewerBinarySetSkipHeader(binaryViewer, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
#if 0
    err = VecView(cellVec, binaryViewer);PYLITH_CHECK_ERROR(err);
#else
    err = PetscObjectTypeCompare((PetscObject) cellVec, VECSEQ, &isseq);PYLITH_CHECK_ERROR(err);
    if (isseq) {err = VecView_Seq(cellVec, binaryViewer);PYLITH_CHECK_ERROR(err);}
    else       {err = VecView_MPI(cellVec, binaryViewer);PYLITH_CHECK_ERROR(err);}
#endif
    err = VecDestroy(&cellVec);PYLITH_CHECK_ERROR(err);
    err = PetscViewerDestroy(&binaryViewer);PYLITH_CHECK_ERROR(err);

    // Create external dataset for cells
    if (!commRank) {
      const hsize_t ndims = 2;
      hsize_t dims[ndims];
      dims[0] = numCells;
      dims[1] = numCorners;
      _h5->createDatasetRawExternal("/topology", "cells", filenameCells.c_str(), dims, ndims, scalartype);
      const int cellDim = mesh.dimension();
      _h5->writeAttribute("/topology/cells", "cell_dim", (void*)&cellDim, H5T_NATIVE_INT);
    } // if
    
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while opening HDF5 file " << _filename << ".\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Unknown error while opening HDF5 file " << _filename << ".";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_END;
} // open

// ----------------------------------------------------------------------
// Close output files.
void
pylith::meshio::DataWriterHDF5Ext::close(void)
{ // close
  PYLITH_METHOD_BEGIN;

  DataWriter::_context = "";

  if (_h5->isOpen()) {
    _h5->close();
  } // if
  _tstampIndex = 0;
  deallocate();

  int commRank = 0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &commRank);
  if (!commRank) {
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
pylith::meshio::DataWriterHDF5Ext::writeVertexField(const PylithScalar t,
						    topology::Field& field,
						    const topology::Mesh& mesh)
{ // writeVertexField
  PYLITH_METHOD_BEGIN;
  
  assert(_h5);

  try {
    const char* context = DataWriter::_context.c_str();

    PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
    MPI_Comm comm;
    PetscMPIInt commRank;
    PetscErrorCode err;

    err = PetscObjectGetComm((PetscObject) dmMesh, &comm);PYLITH_CHECK_ERROR(err);
    err = MPI_Comm_rank(comm, &commRank);PYLITH_CHECK_ERROR(err);
    field.createScatterWithBC(mesh, "", 0, context);
    field.scatterLocalToGlobal(context);

    PetscViewer binaryViewer;

    const hid_t scalartype = (sizeof(double) == sizeof(PylithScalar)) ? H5T_IEEE_F64BE : H5T_IEEE_F32BE;

    // Create external dataset if necessary
    bool createdExternalDataset = false;
    if (_datasets.find(field.label()) != _datasets.end()) {
      binaryViewer = _datasets[field.label()].viewer;
    } else {
      err = PetscViewerBinaryOpen(comm, _datasetFilename(field.label()).c_str(), FILE_MODE_WRITE, &binaryViewer);PYLITH_CHECK_ERROR(err);
      err = PetscViewerBinarySetSkipHeader(binaryViewer, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
      ExternalDataset dataset;
      dataset.numTimeSteps = 0;
      dataset.viewer = binaryViewer;
      _datasets[field.label()] = dataset;
      
      createdExternalDataset = true;
    } // else
    assert(binaryViewer);

    PetscVec vector = field.vector(context);assert(vector);
#if 0
    err = VecView(vector, binaryViewer);PYLITH_CHECK_ERROR(err);
#else
    PetscBool isseq;
    err = PetscObjectTypeCompare((PetscObject) vector, VECSEQ, &isseq);PYLITH_CHECK_ERROR(err);
    if (isseq) {err = VecView_Seq(vector, binaryViewer);PYLITH_CHECK_ERROR(err);}
    else       {err = VecView_MPI(vector, binaryViewer);PYLITH_CHECK_ERROR(err);}
#endif

    ExternalDataset& datasetInfo = _datasets[field.label()];
    ++datasetInfo.numTimeSteps;

    // Update time stamp in "/time, if necessary.      
    if (!commRank) {
      if (_tstampIndex+1 == datasetInfo.numTimeSteps)
        _writeTimeStamp(t);
    } // if

    // Add dataset to HDF5 file, if necessary
    if (createdExternalDataset) {
      PetscDM dm = NULL;
      PetscSection section = NULL;
      PetscInt dof = 0, n, numVerticesLocal = 0, numVertices, vStart;
      PetscIS globalVertexNumbers = NULL;

      err = VecGetDM(vector, &dm);PYLITH_CHECK_ERROR(err);assert(dm);
      err = DMGetDefaultSection(dm, &section);PYLITH_CHECK_ERROR(err);assert(section);

      err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, NULL);PYLITH_CHECK_ERROR(err);
      err = DMPlexGetVertexNumbering(dmMesh, &globalVertexNumbers);PYLITH_CHECK_ERROR(err);
      err = ISGetLocalSize(globalVertexNumbers, &n);PYLITH_CHECK_ERROR(err);
      if (n > 0) {
	const PetscInt *indices = NULL;
	err = ISGetIndices(globalVertexNumbers, &indices);PYLITH_CHECK_ERROR(err);
	err = PetscSectionGetDof(section, vStart, &dof);PYLITH_CHECK_ERROR(err);
	for(PetscInt v = 0; v < n; ++v) {
	  if (indices[v] >= 0) ++numVerticesLocal;
	} // for
	err = ISRestoreIndices(globalVertexNumbers, &indices);PYLITH_CHECK_ERROR(err);
      } // if
      int fiberDimLocal = dof;
      int fiberDim = 0;
      err = MPI_Allreduce(&fiberDimLocal, &fiberDim, 1, MPI_INT, MPI_MAX, comm);PYLITH_CHECK_ERROR(err);
      err = MPI_Allreduce(&numVerticesLocal, &numVertices, 1, MPI_INT, MPI_SUM, comm);PYLITH_CHECK_ERROR(err);
      assert(fiberDim > 0);assert(numVertices > 0);

      datasetInfo.numPoints = numVertices;
      datasetInfo.fiberDim = fiberDim;

      if (!commRank) {
        // Add new external dataset to HDF5 file.
        const int numTimeSteps = DataWriter::_numTimeSteps;
	const hsize_t ndims = 3;
        hsize_t maxDims[ndims];
	if (numTimeSteps > 0) {
          maxDims[0] = H5S_UNLIMITED;
          maxDims[1] = datasetInfo.numPoints;
          maxDims[2] = datasetInfo.fiberDim;
        } else {
          maxDims[0] = 1;
          maxDims[1] = datasetInfo.numPoints;
          maxDims[2] = datasetInfo.fiberDim;
        } // else
        // Create 'vertex_fields' group if necessary.
        if (!_h5->hasGroup("/vertex_fields"))
          _h5->createGroup("/vertex_fields");
	
        _h5->createDatasetRawExternal("/vertex_fields", field.label(), _datasetFilename(field.label()).c_str(), maxDims, ndims, scalartype);
        std::string fullName = std::string("/vertex_fields/") + field.label();
	const char* sattr = topology::FieldBase::vectorFieldString(field.vectorFieldType());
        _h5->writeAttribute(fullName.c_str(), "vector_field_type", sattr);
      } // if
    } else if (!commRank) {
        // Update number of time steps in external dataset info in HDF5 file.
        const int totalNumTimeSteps = DataWriter::_numTimeSteps;
        assert(totalNumTimeSteps > 0);
	
        const hsize_t ndims = 3;
        hsize_t dims[3];
        dims[0] = datasetInfo.numTimeSteps; // update to current value
        dims[1] = datasetInfo.numPoints;
        dims[2] = datasetInfo.fiberDim;
        _h5->extendDatasetRawExternal("/vertex_fields", field.label(), dims, ndims);
    } // if/else

  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " for HDF5 file '" << _filename << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " for HDF5 file '" << _filename << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_END;
} // writeVertexField

// ----------------------------------------------------------------------
// Write field over cells to file.
void
pylith::meshio::DataWriterHDF5Ext::writeCellField(const PylithScalar t,
						  topology::Field& field,
						  const char* label,
						  const int labelId)
{ // writeCellField
  PYLITH_METHOD_BEGIN;

  assert(_h5);

  try {
    const char* context = DataWriter::_context.c_str();

    PetscDM dmMesh = field.mesh().dmMesh();assert(dmMesh);
    MPI_Comm comm;
    PetscMPIInt commRank;
    PetscErrorCode err;

    err = PetscObjectGetComm((PetscObject) dmMesh, &comm);PYLITH_CHECK_ERROR(err);
    err = MPI_Comm_rank(comm, &commRank);PYLITH_CHECK_ERROR(err);
    field.createScatterWithBC(field.mesh(), label ? label : "", labelId, context);
    field.scatterLocalToGlobal(context);

    PetscViewer binaryViewer;

    const hid_t scalartype = (sizeof(double) == sizeof(PylithScalar)) ? H5T_IEEE_F64BE : H5T_IEEE_F32BE;

    // Create external dataset if necessary
    bool createdExternalDataset = false;
    if (_datasets.find(field.label()) != _datasets.end()) {
      binaryViewer = _datasets[field.label()].viewer;
    } else {
      err = PetscViewerBinaryOpen(comm, _datasetFilename(field.label()).c_str(), FILE_MODE_WRITE, &binaryViewer);PYLITH_CHECK_ERROR(err);
      err = PetscViewerBinarySetSkipHeader(binaryViewer, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
      ExternalDataset dataset;
      dataset.numTimeSteps = 0;
      dataset.viewer = binaryViewer;
      _datasets[field.label()] = dataset;

      createdExternalDataset = true;
    } // else
    assert(binaryViewer);

    PetscVec vector = field.vector(context);assert(vector);
#if 0
    err = VecView(vector, binaryViewer);PYLITH_CHECK_ERROR(err);
#else
    PetscBool isseq;
    err = PetscObjectTypeCompare((PetscObject) vector, VECSEQ, &isseq);PYLITH_CHECK_ERROR(err);
    if (isseq) {err = VecView_Seq(vector, binaryViewer);PYLITH_CHECK_ERROR(err);}
    else       {err = VecView_MPI(vector, binaryViewer);PYLITH_CHECK_ERROR(err);}
#endif

    ExternalDataset& datasetInfo = _datasets[field.label()];
    ++datasetInfo.numTimeSteps;

    // Update time stamp in "/time, if necessary.      
    if (!commRank) {
      if (_tstampIndex+1 == datasetInfo.numTimeSteps)
        _writeTimeStamp(t);
    } // if

    // Add dataset to HDF5 file, if necessary
    if (createdExternalDataset) {
      // Get cell information
      PetscSection section = field.localSection();assert(section);
      PetscInt dof = 0, numLocalCells = 0, numCells, cellHeight, cStart, cEnd;
      PetscIS globalCellNumbers;
    
      err = DMPlexGetVTKCellHeight(dmMesh, &cellHeight);PYLITH_CHECK_ERROR(err);
      err = DMPlexGetHeightStratum(dmMesh, cellHeight, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
      if (label) {
	topology::StratumIS cellsIS(dmMesh, label, labelId);
	numCells = cellsIS.size();
	const PetscInt* cells = (numCells > 0) ? cellsIS.points() : 0;
	for(PetscInt c = 0; c < numCells; ++c) {
	  err = PetscSectionGetDof(section, cells[0], &dof);PYLITH_CHECK_ERROR(err);
	  if ((cells[c] >= cStart) && (cells[c] < cEnd)) {
	    ++numLocalCells;
	  } // if
	} // for
      } else {
	err = DMPlexGetCellNumbering(dmMesh, &globalCellNumbers);PYLITH_CHECK_ERROR(err);
	err = ISGetLocalSize(globalCellNumbers, &numCells);PYLITH_CHECK_ERROR(err);
	if (numCells > 0) {
	  const PetscInt *indices = NULL;
	  err = ISGetIndices(globalCellNumbers, &indices);PYLITH_CHECK_ERROR(err);
	  err = PetscSectionGetDof(section, cStart, &dof);PYLITH_CHECK_ERROR(err);
	  for(PetscInt v = 0; v < numCells; ++v) {
	    if (indices[v] >= 0) ++numLocalCells;
	  } // for
	  err = ISRestoreIndices(globalCellNumbers, &indices);PYLITH_CHECK_ERROR(err);
	} // if
      } // if/else
      int fiberDimLocal = dof;
      int fiberDim = 0;
      MPI_Allreduce(&fiberDimLocal, &fiberDim, 1, MPI_INT, MPI_MAX, comm);
      err = MPI_Allreduce(&numLocalCells, &numCells, 1, MPI_INT, MPI_SUM, comm);PYLITH_CHECK_ERROR(err);
      assert(fiberDim > 0);assert(numCells > 0);

      datasetInfo.numPoints = numCells;
      datasetInfo.fiberDim = fiberDim;

      if (!commRank) {
	// Add new external dataset to HDF5 file.	
        const int numTimeSteps = DataWriter::_numTimeSteps;
	const hsize_t ndims = 3;
        hsize_t maxDims[ndims];
        if (numTimeSteps > 0) {
          maxDims[0] = H5S_UNLIMITED;
          maxDims[1] = datasetInfo.numPoints;
          maxDims[2] = datasetInfo.fiberDim;
        } else {
          maxDims[0] = 1;
          maxDims[1] = datasetInfo.numPoints;
          maxDims[2] = fiberDim;
        } // else
        // Create 'cell_fields' group if necessary.
        if (!_h5->hasGroup("/cell_fields"))
          _h5->createGroup("/cell_fields");
	
        _h5->createDatasetRawExternal("/cell_fields", field.label(), _datasetFilename(field.label()).c_str(), maxDims, ndims, scalartype);
        std::string fullName = std::string("/cell_fields/") + field.label();
	const char* sattr = topology::FieldBase::vectorFieldString(field.vectorFieldType());
        _h5->writeAttribute(fullName.c_str(), "vector_field_type", sattr);
      } // if

    } else if (!commRank) {
      // Update number of time steps in external dataset info in HDF5 file.
      const int totalNumTimeSteps = DataWriter::_numTimeSteps;assert(totalNumTimeSteps > 0);
	
      const hsize_t ndims = 3;
      hsize_t dims[3];
      dims[0] = datasetInfo.numTimeSteps; // update to current value
      dims[1] = datasetInfo.numPoints;
      dims[2] = datasetInfo.fiberDim;
      _h5->extendDatasetRawExternal("/cell_fields", field.label(), dims, ndims);
    } // if/else
      
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " for HDF5 file '" << _filename << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " for HDF5 file '" << _filename << "'.";
    throw std::runtime_error(msg.str());
  } // try/catch

  PYLITH_METHOD_END;
} // writeCellField

// ----------------------------------------------------------------------
// Write dataset with names of points to file.
void
pylith::meshio::DataWriterHDF5Ext::writePointNames(const char* const* names,
						   const int numNames,
						   const topology::Mesh& mesh)
{ // writePointNames
  PYLITH_METHOD_BEGIN;

  assert(_h5);

  try {
    PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
    MPI_Comm comm;
    PetscMPIInt commRank;
    PetscErrorCode err = PetscObjectGetComm((PetscObject) dmMesh, &comm);PYLITH_CHECK_ERROR(err);
    err = MPI_Comm_rank(comm, &commRank);PYLITH_CHECK_ERROR(err);

    if (!commRank) {
      _h5->writeDataset("/", "stations", names, numNames);
    } // if
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
pylith::meshio::DataWriterHDF5Ext::_hdf5Filename(void) const
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
// Generate filename for external dataset file.
std::string
pylith::meshio::DataWriterHDF5Ext::_datasetFilename(const char* field) const
{ // _datasetFilename
  PYLITH_METHOD_BEGIN;

  std::ostringstream filenameS;
  std::string filenameH5 = _hdf5Filename();
  const int indexExt = filenameH5.find(".h5");
  filenameS << std::string(filenameH5, 0, indexExt) << "_" << field << ".dat";

  PYLITH_METHOD_RETURN(std::string(filenameS.str()));
} // _datasetFilename

// ----------------------------------------------------------------------
// Write time stamp to file.
void
pylith::meshio::DataWriterHDF5Ext::_writeTimeStamp(const PylithScalar t)
{ // _writeTimeStamp
  PYLITH_METHOD_BEGIN;

  assert(_h5);

  const int ndims = 3;
  const hid_t scalartype = (sizeof(double) == sizeof(PylithScalar)) ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;

  // Each time stamp has a size of 1.
  hsize_t dimsChunk[3]; // Use 3 dims for compatibility with PETSc viewer
  dimsChunk[0] = 1;
  dimsChunk[1] = 1;
  dimsChunk[2] = 1;

  if (!_h5->hasDataset("/time")) {
    // Create dataset
    // Dataset has unknown size.
    hsize_t dims[3];
    dims[0] = H5S_UNLIMITED;
    dims[1] = 1;
    dims[2] = 1;
    _h5->createDataset("/", "time", dims, dimsChunk, ndims, scalartype);
  } // if
  
  // Write time stamp as chunk to HDF5 file.
  // Current dimensions of dataset.
  hsize_t dims[3];
  dims[0] = _tstampIndex+1;
  dims[1] = 1;
  dims[2] = 1;
  const PylithScalar tDim = t * DataWriter::_timeScale;
  _h5->writeDatasetChunk("/", "time", &tDim, dims, dimsChunk, ndims, _tstampIndex, scalartype);
  
  _tstampIndex++;

  PYLITH_METHOD_END;
} // _writeTimeStamp


// End of file 
