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

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::DataWriterHDF5Ext(void) :
  _filename("output.h5"),
  _h5(new HDF5),
  _tstampIndex(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::~DataWriterHDF5Ext(void)
{ // destructor
  delete _h5; _h5 = 0;
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::deallocate(void)
{ // deallocate
  DataWriter<mesh_type, field_type>::deallocate();

  PetscErrorCode err = 0;
  const typename dataset_type::const_iterator& dEnd = _datasets.end();
  for (typename dataset_type::iterator d_iter=_datasets.begin();
       d_iter != dEnd;
       ++d_iter) {
    err = PetscViewerDestroy(&d_iter->second.viewer);CHECK_PETSC_ERROR(err);
  } // for
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor.
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::DataWriterHDF5Ext(const DataWriterHDF5Ext<mesh_type, field_type>& w) :
  DataWriter<mesh_type, field_type>(w),
  _filename(w._filename),
  _h5(new HDF5),
  _tstampIndex(0)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Prepare for writing files.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::open(
						const mesh_type& mesh,
						const int numTimeSteps,
						const char* label,
						const int labelId)
{ // open
  typedef typename mesh_type::SieveMesh SieveMesh;
  typedef typename mesh_type::SieveMesh::label_sequence label_sequence;
  typedef typename mesh_type::SieveMesh::numbering_type numbering_type;
  typedef typename mesh_type::SieveMesh::sieve_type sieve_type;

  assert(_h5);
  _datasets.clear();

  try {
    DataWriter<mesh_type, field_type>::open(mesh, numTimeSteps, label, labelId);
    const char* context = DataWriter<mesh_type, field_type>::_context.c_str();

    DM dmMesh = mesh.dmMesh();
    assert(dmMesh);
    PetscMPIInt    commRank;
    PetscErrorCode err = MPI_Comm_rank(((PetscObject) dmMesh)->comm, &commRank);CHECK_PETSC_ERROR(err);

    if (!commRank) {
      _h5->open(_hdf5Filename().c_str(), H5F_ACC_TRUNC);

      // Create groups
      _h5->createGroup("/topology");
      _h5->createGroup("/geometry");
    } // if
    _tstampIndex = 0;

    PetscViewer binaryViewer;
    
    const hid_t scalartype = (sizeof(double) == sizeof(PylithScalar)) ? 
      H5T_IEEE_F64BE : H5T_IEEE_F32BE;

    // Write vertex coordinates
    const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
    assert(cs);
#if 0
    const ALE::Obj<typename mesh_type::RealSection>& coordinatesSection = 
      sieveMesh->hasRealSection("coordinates_dimensioned") ?
      sieveMesh->getRealSection("coordinates_dimensioned") :
      sieveMesh->getRealSection("coordinates");
    assert(!coordinatesSection.isNull());
    topology::FieldBase::Metadata metadata;
    // :KLUDGE: We would like to use field_type for the coordinates
    // field. However, the mesh coordinates are Field<mesh_type> and
    // field_type can be Field<Mesh> (e.g., displacement field over a
    // SubMesh).
    topology::Field<mesh_type> coordinates(mesh, coordinatesSection, metadata);
    coordinates.label("vertices");
    ALE::Obj<numbering_type> vNumbering = 
      sieveMesh->hasLabel("censored depth") ?
      sieveMesh->getFactory()->getNumbering(sieveMesh, "censored depth", 0) :
      sieveMesh->getFactory()->getNumbering(sieveMesh, 0);
    assert(!vNumbering.isNull());
    coordinates.createScatterWithBC(mesh, vNumbering, context);
    coordinates.scatterSectionToVector(context);
    PetscVec coordinatesVector = coordinates.vector(context);
    assert(coordinatesVector);
#else
    PetscSection coordSection;
    Vec          coordinates, coordinatesVector;
    PetscReal    lengthScale;
    PetscInt     vStart, vEnd, vMax, verticesSize, globalSize, dim, dimLocal = 0;

    /* TODO Get rid of this and use the createScatterWithBC(numbering) code */
    err = DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale);CHECK_PETSC_ERROR(err);
    err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
    err = DMGetCoordinatesLocal(dmMesh, &coordinates);CHECK_PETSC_ERROR(err);
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
    err = DMPlexGetVTKBounds(dmMesh, PETSC_NULL, &vMax);CHECK_PETSC_ERROR(err);
    if (vMax >= 0) {vEnd = PetscMin(vEnd, vMax);}
    for(PetscInt vertex = vStart; vertex < vEnd; ++vertex) {
      err = PetscSectionGetDof(coordSection, vertex, &dimLocal);CHECK_PETSC_ERROR(err);
      if (dimLocal) break;
    }
    err = MPI_Allreduce(&dimLocal, &dim, 1, MPIU_INT, MPI_MAX, mesh.comm());CHECK_PETSC_ERROR(err);
    verticesSize = vEnd - vStart;

    PetscVec     coordVec;
    PetscScalar *coords, *c;

    err = VecCreate(mesh.comm(), &coordVec);CHECK_PETSC_ERROR(err);
    err = VecSetSizes(coordVec, verticesSize*dim, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(coordVec, dim);CHECK_PETSC_ERROR(err);
    err = VecSetFromOptions(coordVec);CHECK_PETSC_ERROR(err);
    err = VecGetArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
    err = VecGetArray(coordinates, &c);CHECK_PETSC_ERROR(err);
    for(PetscInt v = 0; v < vEnd - vStart; ++v) {
      for(PetscInt d = 0; d < dim; ++d) {
          coords[v*dim+d] = c[v*dim+d];
      }
    }
    err = VecRestoreArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
    err = VecRestoreArray(coordinates, &c);CHECK_PETSC_ERROR(err);
    err = VecScale(coordVec, lengthScale);CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) coordVec, "vertices");CHECK_PETSC_ERROR(err);
    err = VecGetSize(coordVec, &globalSize);CHECK_PETSC_ERROR(err);
    globalSize /= cs->spaceDim();
    coordinatesVector = coordVec;
#endif

    const std::string& filenameVertices = _datasetFilename("vertices");
    err = PetscViewerBinaryOpen(((PetscObject) dmMesh)->comm, filenameVertices.c_str(),
				FILE_MODE_WRITE,
				&binaryViewer);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerBinarySetSkipHeader(binaryViewer, PETSC_TRUE);
    CHECK_PETSC_ERROR(err);
    err = VecView(coordinatesVector, binaryViewer); CHECK_PETSC_ERROR(err);
    err = PetscViewerDestroy(&binaryViewer); CHECK_PETSC_ERROR(err);

#if 1
    err = VecDestroy(&coordinatesVector);CHECK_PETSC_ERROR(err);
#endif
    
    // Create external dataset for coordinates    
    if (!commRank) {
      const hsize_t ndims = 2;
      hsize_t dims[ndims];
      dims[0] = globalSize;
      dims[1] = cs->spaceDim();
      _h5->createDatasetRawExternal("/geometry", "vertices", 
				    filenameVertices.c_str(),
				    dims, ndims, scalartype);
    } // if
    
    // Write cells

    // Account for censored cells
#if 0
    int cellDepthLocal = (sieveMesh->depth() == -1) ? -1 : 1;
    int cellDepth = 0;
    err = MPI_Allreduce(&cellDepthLocal, &cellDepth, 1, MPI_INT, MPI_MAX, 
			((PetscObject) dmMesh)->comm);CHECK_PETSC_ERROR(err);
    const int depth = (0 == label) ? cellDepth : labelId;
    const std::string labelName = (0 == label) ?
      ((sieveMesh->hasLabel("censored depth")) ?
       "censored depth" : "depth") : label;
    assert(!sieveMesh->getFactory().isNull());
    const ALE::Obj<numbering_type>& cNumbering = 
      sieveMesh->getFactory()->getNumbering(sieveMesh, labelName, depth);
    assert(!cNumbering.isNull());
    const ALE::Obj<label_sequence>& cells =
      sieveMesh->getLabelStratum(labelName, depth);
    assert(!cells.isNull());
    int numCornersLocal = 0;
    if (cells->size() > 0)
      numCornersLocal = sieveMesh->getNumCellCorners(*cells->begin());
    int numCorners = numCornersLocal;
    err = MPI_Allreduce(&numCornersLocal, &numCorners, 1, MPI_INT, MPI_MAX,
		     ((PetscObject) dmMesh)->comm); CHECK_PETSC_ERROR(err);

    PylithScalar* tmpVertices = 0;
    const int ncells = cNumbering->getLocalSize();
    const int conesSize = ncells*numCorners;
    err = PetscMalloc(sizeof(PylithScalar)*conesSize, &tmpVertices);
    CHECK_PETSC_ERROR(err);

    const Obj<sieve_type>& sieve = sieveMesh->getSieve();
    assert(!sieve.isNull());
    const int closureSize = 
      int(pow(sieve->getMaxConeSize(), sieveMesh->depth()));
    assert(closureSize >= 0);
    ALE::ISieveVisitor::NConeRetriever<sieve_type> ncV(*sieve, closureSize);
    
    int k = 0;
    const typename label_sequence::const_iterator cellsEnd = cells->end();
    for (typename label_sequence::iterator c_iter=cells->begin();
         c_iter != cellsEnd;
         ++c_iter) {
      if (cNumbering->isLocal(*c_iter)) {
        ncV.clear();
        ALE::ISieveTraversal<sieve_type>::orientedClosure(*sieve, *c_iter, ncV);
        const typename ALE::ISieveVisitor::NConeRetriever<sieve_type>::oriented_point_type* cone =
          ncV.getOrientedPoints();
        const int coneSize = ncV.getOrientedSize();
        if (coneSize != numCorners) {
          const char *name;
          std::ostringstream msg;
          err = PetscObjectGetName((PetscObject) dmMesh, &name);CHECK_PETSC_ERROR(err);
          msg << "Inconsistency in topology found for mesh '" << name << "' during output.\n"
              << "Number of vertices (" << coneSize << ") in cell '" << c
              << "' does not expected number of vertices (" << numCorners << ").";
          throw std::runtime_error(msg.str());
        } // if
        for(PetscInt c=0; c < coneSize; ++c)
          tmpVertices[k++] = vNumbering->getIndex(cone[c].first);
      } // if
    }
    PetscVec elemVec;

    err = VecCreateMPIWithArray(((PetscObject) dmMesh)->comm, numCorners, conesSize, PETSC_DETERMINE,
                                tmpVertices, &elemVec); CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) elemVec, "cells");CHECK_PETSC_ERROR(err);
#else
    PetscInt cStart, cEnd, cMax, dof, conesSize, numCells, numCorners, numCornersLocal = 0;

    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
    err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
    err = DMPlexGetVTKBounds(dmMesh, &cMax, PETSC_NULL);CHECK_PETSC_ERROR(err);
    if (cMax >= 0) {cEnd = PetscMin(cEnd, cMax);}
    for(PetscInt cell = cStart; cell < cEnd; ++cell) {
      PetscInt *closure = PETSC_NULL;
      PetscInt  closureSize, v;

      err = DMPlexGetTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
      numCornersLocal = 0;
      for (v = 0; v < closureSize*2; v += 2) {
        if ((closure[v] >= vStart) && (closure[v] < vEnd)) {
          ++numCornersLocal;
        }
      }
      err = DMPlexRestoreTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
      if (numCornersLocal) break;
    }
    err = MPI_Allreduce(&numCornersLocal, &numCorners, 1, MPIU_INT, MPI_MAX, mesh.comm());CHECK_PETSC_ERROR(err);
    if (label) {
      conesSize = 0;
      for(PetscInt cell = cStart; cell < cEnd; ++cell) {
        PetscInt value;

        err = DMPlexGetLabelValue(dmMesh, label, cell, &value);CHECK_PETSC_ERROR(err);
        if (value == labelId) ++conesSize;
      }
      conesSize *= numCorners;
    } else {
      conesSize = (cEnd - cStart)*numCorners;
    }
    CHKMEMA;

    IS              subpointMap, globalVertexNumbers;
    const PetscInt *gvertex = PETSC_NULL, *gpoints = PETSC_NULL;
    PetscVec        cellVec, elemVec;
    PetscScalar    *vertices;

    err = DMPlexGetSubpointMap(dmMesh, &subpointMap);CHECK_PETSC_ERROR(err);
    err = DMPlexGetVertexNumbering(dmMesh, &globalVertexNumbers);CHECK_PETSC_ERROR(err);
    if (subpointMap) {err = ISGetIndices(subpointMap, &gpoints);CHECK_PETSC_ERROR(err);}
    err = ISGetIndices(globalVertexNumbers, &gvertex);CHECK_PETSC_ERROR(err);
    err = VecCreate(mesh.comm(), &cellVec);CHECK_PETSC_ERROR(err);
    err = VecSetSizes(cellVec, conesSize, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(cellVec, numCorners);CHECK_PETSC_ERROR(err);
    err = VecSetFromOptions(cellVec);CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) cellVec, "cells");CHECK_PETSC_ERROR(err);
    err = VecGetArray(cellVec, &vertices);CHECK_PETSC_ERROR(err);
    for(PetscInt cell = cStart, v = 0; cell < cEnd; ++cell) {
      PetscInt *closure = PETSC_NULL;
      PetscInt  closureSize, p;

      if (label) {
        PetscInt value;

        err = DMPlexGetLabelValue(dmMesh, label, cell, &value);CHECK_PETSC_ERROR(err);
        if (value != labelId) continue;
      }
      err = DMPlexGetTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
      for(p = 0; p < closureSize*2; p += 2) {
        if ((closure[p] >= vStart) && (closure[p] < vEnd)) {
          //const PetscInt gv = gpoints ? gpoints[closure[p]] : gvertex[closure[p] - vStart];
          const PetscInt gv = gvertex[closure[p] - vStart];
          vertices[v++] = gv < 0 ? -(gv+1) : gv;
        }
      }
      err = DMPlexRestoreTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
      //assert(v == (cell-cStart+1)*numCorners);
    }
    CHKMEMA;
    err = VecRestoreArray(cellVec, &vertices);CHECK_PETSC_ERROR(err);
    err = VecGetSize(cellVec, &numCells);CHECK_PETSC_ERROR(err);
    numCells /= numCorners;
    elemVec = cellVec;
#endif

    const std::string& filenameCells = _datasetFilename("cells");
    err = PetscViewerBinaryOpen(((PetscObject) dmMesh)->comm, filenameCells.c_str(),
                                FILE_MODE_WRITE, &binaryViewer);CHECK_PETSC_ERROR(err);
    err = PetscViewerBinarySetSkipHeader(binaryViewer, PETSC_TRUE);CHECK_PETSC_ERROR(err);
    err = VecView(elemVec, binaryViewer);CHECK_PETSC_ERROR(err);
    err = VecDestroy(&elemVec);CHECK_PETSC_ERROR(err);
#if 0
    err = PetscFree(tmpVertices);CHECK_PETSC_ERROR(err);
#endif
    err = PetscViewerDestroy(&binaryViewer);CHECK_PETSC_ERROR(err);

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
    msg << "Unknown error while opening HDF5 file " << _filename << ".\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // open

// ----------------------------------------------------------------------
// Close output files.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::close(void)
{ // close
  DataWriter<mesh_type, field_type>::_context = "";

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
} // close

// ----------------------------------------------------------------------
// Write field over vertices to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::writeVertexField(
				            const PylithScalar t,
					    field_type& field,
					    const mesh_type& mesh)
{ // writeVertexField
  typedef typename mesh_type::SieveMesh::numbering_type numbering_type;

  assert(_h5);

  try {
    const char* context = DataWriter<mesh_type, field_type>::_context.c_str();

    DM             dmMesh = mesh.dmMesh();
    PetscMPIInt    commRank;
    PetscErrorCode err;

    assert(dmMesh);
    err = MPI_Comm_rank(((PetscObject) dmMesh)->comm, &commRank);CHECK_PETSC_ERROR(err);
    field.createScatterWithBC(mesh, "", 0, context);
    field.scatterSectionToVector(context);
    PetscVec vector = field.vector(context);
    assert(vector);

    PetscViewer binaryViewer;

    const hid_t scalartype = (sizeof(double) == sizeof(PylithScalar)) ? 
      H5T_IEEE_F64BE : H5T_IEEE_F32BE;

    // Create external dataset if necessary
    bool createdExternalDataset = false;
    if (_datasets.find(field.label()) != _datasets.end()) {
      binaryViewer = _datasets[field.label()].viewer;
    } else {
      err = PetscViewerBinaryOpen(((PetscObject) dmMesh)->comm, 
                                  _datasetFilename(field.label()).c_str(),
                                  FILE_MODE_WRITE, &binaryViewer);
      CHECK_PETSC_ERROR(err);
      err = PetscViewerBinarySetSkipHeader(binaryViewer, PETSC_TRUE);CHECK_PETSC_ERROR(err);
      ExternalDataset dataset;
      dataset.numTimeSteps = 0;
      dataset.viewer = binaryViewer;
      _datasets[field.label()] = dataset;
      
      createdExternalDataset = true;
    } // else
    assert(binaryViewer);

    err = VecView(vector, binaryViewer);CHECK_PETSC_ERROR(err);
    ++_datasets[field.label()].numTimeSteps;

    PetscSection section = field.petscSection();
    PetscInt     dof     = 0, n, numLocalVertices = 0, numVertices, vStart;
    IS           globalVertexNumbers;

    assert(section);
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, PETSC_NULL);CHECK_PETSC_ERROR(err);
    err = DMPlexGetVertexNumbering(dmMesh, &globalVertexNumbers);CHECK_PETSC_ERROR(err);
    err = ISGetSize(globalVertexNumbers, &n);CHECK_PETSC_ERROR(err);
    if (n > 0) {
      const PetscInt *indices;
      err = ISGetIndices(globalVertexNumbers, &indices);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetDof(section, vStart, &dof);CHECK_PETSC_ERROR(err);
      for(PetscInt v = 0; v < n; ++v) {
        if (indices[v] >= 0) ++numLocalVertices;
      }
      err = ISRestoreIndices(globalVertexNumbers, &indices);CHECK_PETSC_ERROR(err);
    }
    int fiberDimLocal = dof;
    int fiberDim = 0;
    err = MPI_Allreduce(&fiberDimLocal, &fiberDim, 1, MPI_INT, MPI_MAX, ((PetscObject) dmMesh)->comm);CHECK_PETSC_ERROR(err);
    err = MPI_Allreduce(&numLocalVertices, &numVertices, 1, MPI_INT, MPI_SUM, ((PetscObject) dmMesh)->comm);CHECK_PETSC_ERROR(err);
    assert(fiberDim > 0);assert(numVertices > 0);

    if (!commRank) {
      if (createdExternalDataset) {
        // Add new external dataset to HDF5 file.
        const int numTimeSteps = DataWriter<mesh_type, field_type>::_numTimeSteps;
        const hsize_t ndims = (numTimeSteps > 0) ? 3 : 2;
        hsize_t maxDims[3];
        if (3 == ndims) {
          maxDims[0] = H5S_UNLIMITED;
          maxDims[1] = numVertices;
          maxDims[2] = fiberDim;
        } else {
          maxDims[0] = numVertices;
          maxDims[1] = fiberDim;
        } // else
        // Create 'vertex_fields' group if necessary.
        if (!_h5->hasGroup("/vertex_fields"))
          _h5->createGroup("/vertex_fields");
	
        _h5->createDatasetRawExternal("/vertex_fields", field.label(),
                                      _datasetFilename(field.label()).c_str(),
                                      maxDims, ndims, scalartype);
        std::string fullName = std::string("/vertex_fields/") + field.label();
        _h5->writeAttribute(fullName.c_str(), "vector_field_type",
                            topology::FieldBase::vectorFieldString(field.vectorFieldType()));
      } else {
        // Update number of time steps in external dataset info in HDF5 file.
        const int totalNumTimeSteps = 
          DataWriter<mesh_type, field_type>::_numTimeSteps;
        assert(totalNumTimeSteps > 0);
        const int numTimeSteps = _datasets[field.label()].numTimeSteps;
	
        const hsize_t ndims = 3;
        hsize_t dims[3];
        dims[0] = numTimeSteps; // update to current value
        dims[1] = numVertices;
        dims[2] = fiberDim;
        _h5->extendDatasetRawExternal("/vertex_fields", field.label(),
                                      dims, ndims);
      } // if/else
      // Update time stamp in "/time, if necessary.
      if (_tstampIndex+1 == _datasets[field.label()].numTimeSteps)
        _writeTimeStamp(t);
    } // if

  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " for HDF5 file '" << _filename << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " for HDF5 file '" << _filename << "'.\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // writeVertexField

// ----------------------------------------------------------------------
// Write field over cells to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::writeCellField(
				       const PylithScalar t,
				       field_type& field,
				       const char* label,
				       const int labelId)
{ // writeCellField
  typedef typename mesh_type::SieveMesh::numbering_type numbering_type;

  assert(_h5);

  try {
    const char* context = DataWriter<mesh_type, field_type>::_context.c_str();

    DM             dmMesh = field.mesh().dmMesh();
    PetscMPIInt    commRank;
    PetscErrorCode err;

    assert(dmMesh);
    err = MPI_Comm_rank(((PetscObject) dmMesh)->comm, &commRank);CHECK_PETSC_ERROR(err);
    field.createScatterWithBC(field.mesh(), label ? label : "", labelId, context);
    field.scatterSectionToVector(context);
    PetscVec vector = field.vector(context);
    assert(vector);

    PetscViewer binaryViewer;

    const hid_t scalartype = (sizeof(double) == sizeof(PylithScalar)) ? 
      H5T_IEEE_F64BE : H5T_IEEE_F32BE;

    // Create external dataset if necessary
    bool createdExternalDataset = false;
    if (_datasets.find(field.label()) != _datasets.end()) {
      binaryViewer = _datasets[field.label()].viewer;
    } else {
      err = PetscViewerBinaryOpen(((PetscObject) dmMesh)->comm,
                                  _datasetFilename(field.label()).c_str(),
                                  FILE_MODE_WRITE, &binaryViewer);
      CHECK_PETSC_ERROR(err);
      err = PetscViewerBinarySetSkipHeader(binaryViewer, PETSC_TRUE);CHECK_PETSC_ERROR(err);
      ExternalDataset dataset;
      dataset.numTimeSteps = 0;
      dataset.viewer = binaryViewer;
      _datasets[field.label()] = dataset;

      createdExternalDataset = true;
    } // else
    assert(binaryViewer);

    err = VecView(vector, binaryViewer);CHECK_PETSC_ERROR(err);
    ++_datasets[field.label()].numTimeSteps;

    PetscSection section = field.petscSection();
    PetscInt     dof     = 0, n, numLocalCells = 0, numCells, cellHeight, cStart, cEnd;
    IS           globalCellNumbers;

    assert(section);
    err = DMPlexGetVTKCellHeight(dmMesh, &cellHeight);CHECK_PETSC_ERROR(err);
    err = DMPlexGetHeightStratum(dmMesh, cellHeight, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
    if (label) {
      IS pointIS;

      DMPlexGetStratumIS(dmMesh, label, labelId, &pointIS);CHECK_PETSC_ERROR(err);
      err = ISGetLocalSize(pointIS, &n);CHECK_PETSC_ERROR(err);
      if (n > 0) {
        const PetscInt *indices;
        err = ISGetIndices(pointIS, &indices);CHECK_PETSC_ERROR(err);
        for(PetscInt c = 0; c < n; ++c) {
          if ((indices[c] >= cStart) && (indices[c] < cEnd)) {
            err = PetscSectionGetDof(section, indices[0], &dof);CHECK_PETSC_ERROR(err);
            ++numLocalCells;
          }
        }
        err = ISRestoreIndices(pointIS, &indices);CHECK_PETSC_ERROR(err);
      }
    } else {
      err = DMPlexGetCellNumbering(dmMesh, &globalCellNumbers);CHECK_PETSC_ERROR(err);
      err = ISGetLocalSize(globalCellNumbers, &n);CHECK_PETSC_ERROR(err);
      if (n > 0) {
        const PetscInt *indices;
        err = ISGetIndices(globalCellNumbers, &indices);CHECK_PETSC_ERROR(err);
        err = PetscSectionGetDof(section, cStart, &dof);CHECK_PETSC_ERROR(err);
        for(PetscInt v = 0; v < n; ++v) {
          if (indices[v] >= 0) ++numLocalCells;
        }
        err = ISRestoreIndices(globalCellNumbers, &indices);CHECK_PETSC_ERROR(err);
      }
    }
    int fiberDimLocal = dof;
    int fiberDim = 0;
    MPI_Allreduce(&fiberDimLocal, &fiberDim, 1, MPI_INT, MPI_MAX, ((PetscObject) dmMesh)->comm);
    err = MPI_Allreduce(&numLocalCells, &numCells, 1, MPI_INT, MPI_SUM, ((PetscObject) dmMesh)->comm);CHECK_PETSC_ERROR(err);
    assert(fiberDim > 0);assert(numCells > 0);

    if (!commRank) {
      if (createdExternalDataset) {
      // Add new external dataset to HDF5 file.

        const int numTimeSteps =
          DataWriter<mesh_type, field_type>::_numTimeSteps;
        const hsize_t ndims = (numTimeSteps > 0) ? 3 : 2;
        hsize_t maxDims[3];
        if (3 == ndims) {
          maxDims[0] = H5S_UNLIMITED;
          maxDims[1] = numCells;
          maxDims[2] = fiberDim;
        } else {
          maxDims[0] = numCells;
          maxDims[1] = fiberDim;
        } // else
        // Create 'cell_fields' group if necessary.
        if (!_h5->hasGroup("/cell_fields"))
          _h5->createGroup("/cell_fields");
	
        _h5->createDatasetRawExternal("/cell_fields", field.label(),
                                      _datasetFilename(field.label()).c_str(),
                                      maxDims, ndims, scalartype);
        std::string fullName = std::string("/cell_fields/") + field.label();
        _h5->writeAttribute(fullName.c_str(), "vector_field_type",
                            topology::FieldBase::vectorFieldString(field.vectorFieldType()));
      } else {
        // Update number of time steps in external dataset info in HDF5 file.
        const int totalNumTimeSteps = 
          DataWriter<mesh_type, field_type>::_numTimeSteps;
        assert(totalNumTimeSteps > 0);
        const int numTimeSteps = _datasets[field.label()].numTimeSteps;
	
        const hsize_t ndims = 3;
        hsize_t dims[3];
        dims[0] = numTimeSteps; // update to current value
        dims[1] = numCells;
        dims[2] = fiberDim;
        _h5->extendDatasetRawExternal("/cell_fields", field.label(), dims, ndims);
      } // if/else
      // Update time stamp in "/time, if necessary.
      if (_tstampIndex+1 == _datasets[field.label()].numTimeSteps)
        _writeTimeStamp(t);
    } // if

  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " for HDF5 file '" << _filename << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " for HDF5 file '" << _filename << "'.\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // writeCellField

// ----------------------------------------------------------------------
// Generate filename for HDF5 file.
template<typename mesh_type, typename field_type>
std::string
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::_hdf5Filename(void) const
{ // _hdf5Filename
  std::ostringstream filename;
  const int indexExt = _filename.find(".h5");
  const int numTimeSteps = DataWriter<mesh_type, field_type>::_numTimeSteps;
  if (0 == numTimeSteps) {
    filename << std::string(_filename, 0, indexExt) << "_info.h5";
  } else {
    filename << _filename;
  } // if/else

  return std::string(filename.str());
} // _hdf5Filename

// ----------------------------------------------------------------------
// Generate filename for external dataset file.
template<typename mesh_type, typename field_type>
std::string
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::_datasetFilename(const char* field) const
{ // _datasetFilename
  std::ostringstream filenameS;
  std::string filenameH5 = _hdf5Filename();
  const int indexExt = filenameH5.find(".h5");
  filenameS << std::string(filenameH5, 0, indexExt) << "_" << field << ".dat";

  return std::string(filenameS.str());
} // _datasetFilename

// ----------------------------------------------------------------------
// Write time stamp to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::_writeTimeStamp(
						  const PylithScalar t)
{ // _writeTimeStamp
  assert(_h5);

  const int ndims = 3;
  const hid_t scalartype = (sizeof(double) == sizeof(PylithScalar)) ? 
    H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;

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
    _h5->createDataset("/", "time", dims, dimsChunk, ndims, 
		       scalartype);
  } // if
  
  // Write time stamp as chunk to HDF5 file.
  // Current dimensions of dataset.
  hsize_t dims[3];
  dims[0] = _tstampIndex+1;
  dims[1] = 1;
  dims[2] = 1;
  const PylithScalar tDim = t * DataWriter<mesh_type, field_type>::_timeScale;
  _h5->writeDatasetChunk("/", "time", &tDim, dims, dimsChunk, ndims, 
			 _tstampIndex, scalartype);
  
  _tstampIndex++;
} // _writeTimeStamp



// End of file 
