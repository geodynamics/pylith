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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "HDF5.hh" // USES HDF5

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::DataWriterHDF5Ext(void) :
  _filename("output.h5"),
  _h5(new HDF5)
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
       ++d_iter)
    if (d_iter->second.viewer) {
      PetscViewerDestroy(d_iter->second.viewer);
      d_iter->second.viewer = 0;
    } // if
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor.
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::DataWriterHDF5Ext(const DataWriterHDF5Ext<mesh_type, field_type>& w) :
  DataWriter<mesh_type, field_type>(w),
  _filename(w._filename),
  _h5(new HDF5)
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

    PetscErrorCode err = 0;
    
    _h5->open(_hdf5Filename().c_str(), H5F_ACC_TRUNC);

    // Create groups
    _h5->createGroup("/topology");
    _h5->createGroup("/geometry");
    _h5->createGroup("/vertex_fields");
    _h5->createGroup("/cell_fields");

    PetscViewer binaryViewer;
    
    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
    assert(!sieveMesh.isNull());

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
    ALE::Obj<numbering_type> vNumbering;
    if (sieveMesh->hasLabel("censored depth")) { // Remove Lagrange vertices
      vNumbering = sieveMesh->getFactory()->getNumbering(sieveMesh,
							 "censored depth", 0);
      coordinates.createScatter(vNumbering, context);
    } else {
      coordinates.createScatter(context);
    } // if/else
    coordinates.scatterSectionToVector(context);
    PetscVec coordinatesVector = coordinates.vector(context);
    assert(coordinatesVector);

    const std::string& filenameVertices = _datasetFilename("vertices");
    err = PetscViewerBinaryCreate(sieveMesh->comm(), &binaryViewer);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerBinaryOpen(sieveMesh->comm(), filenameVertices.c_str(),
				FILE_MODE_WRITE,
				&binaryViewer);
    CHECK_PETSC_ERROR(err);
    err = VecView(coordinatesVector, binaryViewer); CHECK_PETSC_ERROR(err);
    err = PetscViewerDestroy(binaryViewer); CHECK_PETSC_ERROR(err);
    binaryViewer = 0;
    
    // Create external dataset for coordinates    
    const hsize_t ndims = 2;
    hsize_t dims[ndims];
    dims[0] = vNumbering->getGlobalSize();
    dims[1] = cs->spaceDim();
    _h5->createDatasetRawExternal("/geometry", "vertices", 
				  filenameVertices.c_str(),
				  dims, ndims, H5T_NATIVE_DOUBLE);
    
    // Write cells

    // Account for censored cells
    const int cellDepth = (sieveMesh->depth() == -1) ? -1 : 1;
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
    err = MPI_Reduce(&numCornersLocal, &numCorners, 1, MPI_INT, MPI_MAX, 0, 
		     sieveMesh->comm()); CHECK_PETSC_ERROR(err);

    PetscScalar* tmpVertices = 0;
    const int ncells = cNumbering->getLocalSize();
    const int conesSize = ncells*numCorners;
    err = PetscMalloc(sizeof(PetscScalar)*conesSize, &tmpVertices);
    CHECK_PETSC_ERROR(err);

    const Obj<sieve_type>& sieve = sieveMesh->getSieve();
    assert(!sieve.isNull());
    ALE::ISieveVisitor::NConeRetriever<sieve_type> 
      ncV(*sieve, (size_t) pow((double) sieve->getMaxConeSize(), 
			       std::max(0, sieveMesh->depth())));

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
          for (int c=0; c < coneSize; ++c)
            tmpVertices[k++] = vNumbering->getIndex(cone[c].first);
      } // if

    Vec elemVec;
    err = VecCreateMPIWithArray(sieveMesh->comm(), conesSize, PETSC_DETERMINE,
				tmpVertices, &elemVec); CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) elemVec,
			     "cells");CHECK_PETSC_ERROR(err);

    const std::string& filenameCells = _datasetFilename("cells");
    err = PetscViewerBinaryCreate(sieveMesh->comm(), &binaryViewer);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerBinaryOpen(sieveMesh->comm(), filenameCells.c_str(),
				FILE_MODE_WRITE,
				&binaryViewer);
    CHECK_PETSC_ERROR(err);
    err = VecView(elemVec, binaryViewer); CHECK_PETSC_ERROR(err);
    err = VecDestroy(elemVec); CHECK_PETSC_ERROR(err);
    err = PetscFree(tmpVertices); CHECK_PETSC_ERROR(err);
    err = PetscViewerDestroy(binaryViewer); CHECK_PETSC_ERROR(err);
    binaryViewer = 0;

    // Create external dataset for cells
    dims[0] = cNumbering->getGlobalSize();
    dims[1] = numCorners;
    _h5->createDatasetRawExternal("/topology", "cells", filenameCells.c_str(),
				  dims, ndims, H5T_NATIVE_DOUBLE);
    
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
  deallocate();
} // close

// ----------------------------------------------------------------------
// Write field over vertices to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::writeVertexField(
				            const double t,
					    field_type& field,
					    const mesh_type& mesh)
{ // writeVertexField
  assert(_h5);

  try {
    const char* context = DataWriter<mesh_type, field_type>::_context.c_str();

    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
    assert(!sieveMesh.isNull());
    const std::string labelName = 
      (sieveMesh->hasLabel("censored depth")) ? "censored depth" : "depth";
    const ALE::Obj<typename mesh_type::SieveMesh::numbering_type>& vNumbering =
      sieveMesh->getFactory()->getNumbering(sieveMesh, labelName, 0);
    assert(!vNumbering.isNull());

    if (sieveMesh->hasLabel("censored depth")) { // Remove Lagrange vertices
      field.createScatter(vNumbering, context);
    } else {
      field.createScatter(context);
    } // if/else
    field.scatterSectionToVector(context);
    PetscVec vector = field.vector(context);
    assert(vector);

    PetscViewer binaryViewer;

    // Create external dataset if necessary
    bool createdExternalDataset = false;
    if (_datasets.find(field.label()) != _datasets.end()) {
      binaryViewer = _datasets[field.label()].viewer;
    } else {
      PetscViewerBinaryCreate(sieveMesh->comm(), &binaryViewer);
      PetscViewerBinaryOpen(sieveMesh->comm(), 
			    _datasetFilename(field.label()).c_str(),
			    FILE_MODE_WRITE, &binaryViewer);
      ExternalDataset dataset;
      dataset.numTimeSteps = 0;
      dataset.viewer = binaryViewer;
      _datasets[field.label()] = dataset;
      
      createdExternalDataset = true;
    } // else

    PetscErrorCode err = VecView(vector, binaryViewer);
    CHECK_PETSC_ERROR(err);
    ++_datasets[field.label()].numTimeSteps;

    const ALE::Obj<typename mesh_type::RealSection>& section = 
      field.section();
    assert(!section.isNull());
    assert(!sieveMesh->getLabelStratum(labelName, 0).isNull());
    const int fiberDimLocal = 
      (sieveMesh->getLabelStratum(labelName, 0)->size() > 0) ? 
      section->getFiberDimension(*sieveMesh->getLabelStratum(labelName, 
							     0)->begin()) : 0;
    int fiberDim = 0;
    MPI_Allreduce((void *) &fiberDimLocal, (void *) &fiberDim, 1, 
		  MPI_INT, MPI_MAX, field.mesh().comm());
    assert(fiberDim > 0);

    if (createdExternalDataset) {
      // Add new external dataset to HDF5 file.
      const int numTimeSteps = DataWriter<mesh_type, field_type>::_numTimeSteps;
      const hsize_t ndims = (numTimeSteps > 0) ? 3 : 2;
      hsize_t* dims = (ndims > 0) ? new hsize_t[ndims] : 0;
      if (3 == ndims) {
	dims[0] = 1; // external file only constains 1 time step so far.
	dims[1] = vNumbering->getGlobalSize();
	dims[2] = fiberDim;
      } else {
	dims[0] = vNumbering->getGlobalSize();
	dims[1] = fiberDim;
      } // else
      _h5->createDatasetRawExternal("/vertex_fields", field.label(),
				    _datasetFilename(field.label()).c_str(),
				    dims, ndims, H5T_NATIVE_DOUBLE);
      delete[] dims; dims = 0;
    } else {
      // Update number of time steps in external dataset info in HDF5 file.
      const int totalNumTimeSteps = 
	DataWriter<mesh_type, field_type>::_numTimeSteps;
      assert(totalNumTimeSteps > 0);
      const int numTimeSteps = _datasets[field.label()].numTimeSteps;

      const hsize_t ndims = 3;
      hsize_t dims[ndims];
      dims[0] = numTimeSteps; // update to current value
      dims[1] = vNumbering->getGlobalSize();
      dims[2] = fiberDim;
      _h5->extendDatasetRawExternal("/vertex_fields", field.label(),
				    dims, ndims);
    } // if/else

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
				       const double t,
				       field_type& field,
				       const char* label,
				       const int labelId)
{ // writeCellField
  assert(_h5);

  try {
    // :TODO: Must account for possible presence of 'censored depth'
    // and censor the appropriate vertices.

    const char* context = DataWriter<mesh_type, field_type>::_context.c_str();

    field.createScatter(context);
    field.scatterSectionToVector(context);
    PetscVec vector = field.vector(context);
    assert(vector);

    // :TODO: Need to account for censored cells.

    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = 
      field.mesh().sieveMesh();
    assert(!sieveMesh.isNull());

    PetscViewer binaryViewer;

    // Create external dataset if necessary
    bool createdExternalDataset = false;
    if (_datasets.find(field.label()) != _datasets.end()) {
      binaryViewer = _datasets[field.label()].viewer;
    } else {
      PetscViewerBinaryCreate(sieveMesh->comm(), &binaryViewer);
      PetscViewerBinaryOpen(sieveMesh->comm(),
			    _datasetFilename(field.label()).c_str(),
			    FILE_MODE_WRITE, &binaryViewer);
      ExternalDataset dataset;
      dataset.numTimeSteps = 0;
      dataset.viewer = binaryViewer;
      _datasets[field.label()] = dataset;

      createdExternalDataset = true;
    } // else

    PetscErrorCode err = VecView(vector, binaryViewer);
    CHECK_PETSC_ERROR(err);
    ++_datasets[field.label()].numTimeSteps;

    if (createdExternalDataset) {
      // Add new external dataset to HDF5 file.
      const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = 
	field.mesh().sieveMesh();
      assert(!sieveMesh.isNull());
      const int cellDepth = (sieveMesh->depth() == -1) ? -1 : 1;
      const int depth = (0 == label) ? cellDepth : labelId;
      const std::string labelName = (0 == label) ?
	((sieveMesh->hasLabel("censored depth")) ?
	 "censored depth" : "depth") : label;
      assert(!sieveMesh->getFactory().isNull());
      const ALE::Obj<typename mesh_type::SieveMesh::numbering_type>& numbering = 
	sieveMesh->getFactory()->getNumbering(sieveMesh, labelName, depth);
      assert(!numbering.isNull());
      assert(!sieveMesh->getLabelStratum(labelName, depth).isNull());
      const ALE::Obj<typename mesh_type::RealSection>& section = field.section();
      assert(!section.isNull());
      
      const int fiberDimLocal = 
	(sieveMesh->getLabelStratum(labelName, depth)->size() > 0) ? 
	section->getFiberDimension(*sieveMesh->getLabelStratum(labelName, depth)->begin()) : 0;
      int fiberDim = 0;
      MPI_Allreduce((void *) &fiberDimLocal, (void *) &fiberDim, 1, 
		    MPI_INT, MPI_MAX, field.mesh().comm());
      assert(fiberDim > 0);

      const int numTimeSteps = DataWriter<mesh_type, field_type>::_numTimeSteps;
      const hsize_t ndims = (numTimeSteps > 0) ? 3 : 2;
      hsize_t* dims = (ndims > 0) ? new hsize_t[ndims] : 0;
      if (3 == ndims) {
	dims[0] = 1; // external file only constains 1 time step so far.
	dims[1] = numbering->getGlobalSize();
	dims[2] = fiberDim;
      } else {
	dims[0] = numbering->getGlobalSize();
	dims[1] = fiberDim;
      } // else
      _h5->createDatasetRawExternal("/cell_fields", field.label(),
				    _datasetFilename(field.label()).c_str(),
				    dims, ndims, H5T_NATIVE_DOUBLE);
      delete[] dims; dims = 0;
    } // else

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


// End of file 
