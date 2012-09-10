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
  const typename dataset_type::const_iterator& dEnd = _datasets.end();
  for (typename dataset_type::iterator d_iter=_datasets.begin();
       d_iter != dEnd;
       ++d_iter) {
    PetscErrorCode err = PetscViewerDestroy(&d_iter->second.viewer);
    CHECK_PETSC_ERROR(err);
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

    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
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
    const ALE::Obj<typename mesh_type::RealSection>& coordinatesSection = 
      sieveMesh->hasRealSection("coordinates_dimensioned") ?
      sieveMesh->getRealSection("coordinates_dimensioned") :
      sieveMesh->getRealSection("coordinates");
    assert(!coordinatesSection.isNull());
    const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
    assert(cs);
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

    const std::string& filenameVertices = _datasetFilename("vertices");
    err = PetscViewerBinaryOpen(((PetscObject) dmMesh)->comm, filenameVertices.c_str(),
				FILE_MODE_WRITE,
				&binaryViewer);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerBinarySetSkipHeader(binaryViewer, PETSC_TRUE);
    CHECK_PETSC_ERROR(err);
    err = VecView(coordinatesVector, binaryViewer); CHECK_PETSC_ERROR(err);
    err = PetscViewerDestroy(&binaryViewer); CHECK_PETSC_ERROR(err);
    
    // Create external dataset for coordinates    
    if (!commRank) {
      const hsize_t ndims = 2;
      hsize_t dims[ndims];
      dims[0] = vNumbering->getGlobalSize();
      dims[1] = cs->spaceDim();
      _h5->createDatasetRawExternal("/geometry", "vertices", 
				    filenameVertices.c_str(),
				    dims, ndims, scalartype);
    } // if
    
    // Write cells

    // Account for censored cells
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
	 ++c_iter)
      if (cNumbering->isLocal(*c_iter)) {
	ncV.clear();
	ALE::ISieveTraversal<sieve_type>::orientedClosure(*sieve, *c_iter, ncV);
	const typename ALE::ISieveVisitor::NConeRetriever<sieve_type>::oriented_point_type* cone =
	  ncV.getOrientedPoints();
	const int coneSize = ncV.getOrientedSize();
	if (coneSize != numCorners) {
	  std::ostringstream msg;
	  msg << "Inconsistency in topology found for mesh '"
	      << sieveMesh->getName() << "' during output.\n"
	      << "Number of vertices (" << coneSize << ") in cell '"
	      << *c_iter << "' does not expected number of vertices ("
	      << numCorners << ").";
	  throw std::runtime_error(msg.str());
	} // if
	for (int c=0; c < coneSize; ++c)
	  tmpVertices[k++] = vNumbering->getIndex(cone[c].first);
      } // if

    PetscVec elemVec;
    err = VecCreateMPIWithArray(((PetscObject) dmMesh)->comm, numCorners, conesSize, PETSC_DETERMINE,
				tmpVertices, &elemVec); CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) elemVec,
			     "cells");CHECK_PETSC_ERROR(err);

    const std::string& filenameCells = _datasetFilename("cells");
    err = PetscViewerBinaryOpen(((PetscObject) dmMesh)->comm, filenameCells.c_str(),
				FILE_MODE_WRITE,
				&binaryViewer);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerBinarySetSkipHeader(binaryViewer, PETSC_TRUE);
    CHECK_PETSC_ERROR(err);
    err = VecView(elemVec, binaryViewer); CHECK_PETSC_ERROR(err);
    err = VecDestroy(&elemVec); CHECK_PETSC_ERROR(err);
    err = PetscFree(tmpVertices); CHECK_PETSC_ERROR(err);
    err = PetscViewerDestroy(&binaryViewer); CHECK_PETSC_ERROR(err);

    // Create external dataset for cells
    if (!commRank) {
      const hsize_t ndims = 2;
      hsize_t dims[ndims];
      dims[0] = cNumbering->getGlobalSize();
      dims[1] = numCorners;
      _h5->createDatasetRawExternal("/topology", "cells", filenameCells.c_str(),
				    dims, ndims, scalartype);
      const int cellDim = mesh.dimension();
      _h5->writeAttribute("/topology/cells", "cell_dim", (void*)&cellDim,
			  H5T_NATIVE_INT);
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

    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = field.mesh().sieveMesh();
    assert(!sieveMesh.isNull());
    DM dmMesh = mesh.dmMesh();
    assert(dmMesh);
    PetscMPIInt    commRank;
    PetscErrorCode err = MPI_Comm_rank(((PetscObject) dmMesh)->comm, &commRank);CHECK_PETSC_ERROR(err);

    const std::string labelName = 
      (sieveMesh->hasLabel("censored depth")) ? "censored depth" : "depth";
    ALE::Obj<numbering_type> vNumbering = 
      sieveMesh->hasLabel("censored depth") ?
      sieveMesh->getFactory()->getNumbering(sieveMesh, "censored depth", 0) :
      sieveMesh->getFactory()->getNumbering(sieveMesh, 0);
    assert(!vNumbering.isNull());
    field.createScatterWithBC(mesh, vNumbering, context);
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
      err = PetscViewerBinarySetSkipHeader(binaryViewer, PETSC_TRUE);
      CHECK_PETSC_ERROR(err);
      ExternalDataset dataset;
      dataset.numTimeSteps = 0;
      dataset.viewer = binaryViewer;
      _datasets[field.label()] = dataset;
      
      createdExternalDataset = true;
    } // else
    assert(binaryViewer);

    err = VecView(vector, binaryViewer);
    CHECK_PETSC_ERROR(err);
    ++_datasets[field.label()].numTimeSteps;

    PetscSection section = field.petscSection();
    PetscInt     dof     = 0;
    assert(section);
    assert(!sieveMesh->getLabelStratum(labelName, 0).isNull());
    if (sieveMesh->getLabelStratum(labelName, 0)->size() > 0) {
      err = PetscSectionGetDof(section, *sieveMesh->getLabelStratum(labelName, 0)->begin(), &dof);CHECK_PETSC_ERROR(err);
    }
    int fiberDimLocal = dof;
    int fiberDim = 0;
    MPI_Allreduce(&fiberDimLocal, &fiberDim, 1, MPI_INT, MPI_MAX, ((PetscObject) dmMesh)->comm);
    assert(fiberDim > 0);

    if (!commRank) {
      if (createdExternalDataset) {
        // Add new external dataset to HDF5 file.
        const int numTimeSteps = DataWriter<mesh_type, field_type>::_numTimeSteps;
        const hsize_t ndims = (numTimeSteps > 0) ? 3 : 2;
        hsize_t maxDims[3];
        if (3 == ndims) {
          maxDims[0] = H5S_UNLIMITED;
          maxDims[1] = vNumbering->getGlobalSize();
          maxDims[2] = fiberDim;
        } else {
          maxDims[0] = vNumbering->getGlobalSize();
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
        dims[1] = vNumbering->getGlobalSize();
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

    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = field.mesh().sieveMesh();
    assert(!sieveMesh.isNull());
    DM dmMesh = field.mesh().dmMesh();
    assert(dmMesh);
    PetscMPIInt    commRank;
    PetscErrorCode err = MPI_Comm_rank(((PetscObject) dmMesh)->comm, &commRank);CHECK_PETSC_ERROR(err);

    int cellDepthLocal = (sieveMesh->depth() == -1) ? -1 : 1;
    int cellDepth = 0;
    err = MPI_Allreduce(&cellDepthLocal, &cellDepth, 1, MPI_INT, MPI_MAX, 
			((PetscObject) dmMesh)->comm);CHECK_PETSC_ERROR(err);
    const int depth = (0 == label) ? cellDepth : labelId;
    const std::string labelName = (0 == label) ?
      ((sieveMesh->hasLabel("censored depth")) ?
       "censored depth" : "depth") : label;
    assert(!sieveMesh->getFactory().isNull());
    const ALE::Obj<typename mesh_type::SieveMesh::numbering_type>& numbering = 
      sieveMesh->getFactory()->getNumbering(sieveMesh, labelName, depth);
    assert(!numbering.isNull());
    field.createScatterWithBC(field.mesh(), numbering, context);
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
      err = PetscViewerBinarySetSkipHeader(binaryViewer, PETSC_TRUE);
      CHECK_PETSC_ERROR(err);
      ExternalDataset dataset;
      dataset.numTimeSteps = 0;
      dataset.viewer = binaryViewer;
      _datasets[field.label()] = dataset;

      createdExternalDataset = true;
    } // else
    assert(binaryViewer);

    err = VecView(vector, binaryViewer);
    CHECK_PETSC_ERROR(err);
    ++_datasets[field.label()].numTimeSteps;

    PetscSection section = field.petscSection();
    PetscInt     dof     = 0;
    assert(section);
    assert(!sieveMesh->getLabelStratum(labelName, 0).isNull());
    if (sieveMesh->getLabelStratum(labelName, 0)->size() > 0) {
      err = PetscSectionGetDof(section, *sieveMesh->getLabelStratum(labelName, 0)->begin(), &dof);CHECK_PETSC_ERROR(err);
    }
    int fiberDimLocal = dof;
    int fiberDim = 0;
    MPI_Allreduce(&fiberDimLocal, &fiberDim, 1, MPI_INT, MPI_MAX, ((PetscObject) dmMesh)->comm);
    assert(fiberDim > 0);

    if (!commRank) {
      if (createdExternalDataset) {
      // Add new external dataset to HDF5 file.

	const int numTimeSteps =
	  DataWriter<mesh_type, field_type>::_numTimeSteps;
	const hsize_t ndims = (numTimeSteps > 0) ? 3 : 2;
	hsize_t maxDims[3];
	if (3 == ndims) {
	  maxDims[0] = H5S_UNLIMITED;
	  maxDims[1] = numbering->getGlobalSize();
	  maxDims[2] = fiberDim;
	} else {
	  maxDims[0] = numbering->getGlobalSize();
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
	dims[1] = numbering->getGlobalSize();
	dims[2] = fiberDim;
	_h5->extendDatasetRawExternal("/cell_fields", field.label(),
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
