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

#include <petscmesh_viewers.hh> // USES HDF5Viewer

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
pylith::meshio::DataWriterHDF5Ext<mesh_type, field_type>::deallocate(void)
{ // deallocate
  const dataset_type::const_iterator& dEnd = _datasets.end();
  for (dataset_type::const_iterator d_iter=_datasets.begin();
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
void
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::open(
						const mesh_type& mesh,
						const int numTimeSteps,
						const char* label =0,
						const int labelId =0)
{ // open
  assert(!_h5);
  _datasets->clear();

  try {
    PetscErrorCode err = 0;
    
    _h5->filename(_hdf5Filename(), H5_ACC_TRUNC);

    // Create groups
    _h5->createGroup("/topology");
    _h5->createGroup("/geometry");
    _h5->createGroup("/vertex_fields");
    _h5->createGroup("/cell_fields");

    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
    assert(!sieveMesh.isNull());

    // Write vertex coordinates
    const char* filenameVertices = _datasetFilename("vertices");
    PetscViewer binaryViewer;
    err = PetscViewerBinaryCreate(sieveMesh->comm(), &binaryViewer);
    PETSC_CHECK_ERROR(err);
    err = PetscViewerBinaryOpen(sieveMesh->comm(), 
				_datasetFilename("vertices").c_str(),
				FILE_MODE_WRITE,
				&binaryViewer);

    const ALE::Obj<typename mesh_type::RealSection>& coordinatesSection = 
      sieveMesh->getRealSection("coordinates_dimensioned");

    topology::FieldBase::Metadata metadata;
    // :KLUDGE: We would like to use field_type for the coordinates
    // field. However, the mesh coordinates are Field<mesh_type> and
    // field_type can be Field<Mesh> (e.g., displacement field over a
    // SubMesh).
    topology::Field<mesh_type> coordinates(mesh, coordinatesSection, metadata);
    coordinates.label("vertices");
    // :TODO: Need to test for presence of 'censored depth' label and
    // use it to censor vertices. Need to create Vec consistent with
    // censored vertices. See DataWriterVTK and VTKViewer for how
    // 'censored depth' is used.
    coordinates.createVector();
    coordinates.createScatter();
    coordinates.scatterSectionToVector();
    err = VecView(coordinates.vector(), binaryViewer); CHECK_PETSC_ERROR(err);
    err = PetscViewerDestroy(binaryViewer); CHECK_PETSC_ERROR(err);
    binaryViewer = 0;

    // Create external dataset for coordinates

    // :TODO: Update this to use sizes from numbering to account for
    // censored vertices.
    const ALE::Obj<typename mesh_type::SieveMesh::label_sequence>& vertices =
      sieveMesh->depthStratum(0);
    assert(!vertices.isNull());
    const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
    const hsize_t ndims = 2;
    hsize_t dims[ndims];
    dims[0] = vertices->size();
    dims[1] = cs->spaceDim(();
    _h5->createDatasetRawExternal("/geometry", "vertices", filenameVertices,
				  dims, ndims, H5_NATIVE_DOUBLE);

#if 0
    // Write cells

    // Create external dataset for cells

    Vec          elemVec;
    PetscScalar *tmpVertices;
    PetscBool    columnMajor = PETSC_FALSE;

    typedef ALE::OrientedConeSectionV<typename mesh_type::SieveMesh::sieve_type> oriented_cones_wrapper_type;
    Obj<oriented_cones_wrapper_type> cones = new oriented_cones_wrapper_type(sieveMesh->getSieve());

    // Hack right now, move to HDF5 Section viewer
    err = PetscMalloc(sizeof(PetscScalar)*cones->size(), &tmpVertices);CHECK_PETSC_ERROR(err);
    for(int p = sieveMesh->getSieve()->getChart().min(), i = 0; p < sieveMesh->getSieve()->getChart().max(); ++p) {
      const int coneSize = cones->getFiberDimension(p);
      const typename oriented_cones_wrapper_type::value_type *vertices = cones->restrictPoint(p);

      for(int c = 0; c < coneSize; ++c, ++i) {
        tmpVertices[i] = vertices[c].first;
      }
    }
    err = VecCreateMPIWithArray(sieveMesh->comm(), cones->size(), PETSC_DETERMINE, tmpVertices, &elemVec);CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) elemVec, "cells");CHECK_PETSC_ERROR(err);
    err = VecView(elemVec, _viewer);CHECK_PETSC_ERROR(err);
    err = VecDestroy(elemVec);CHECK_PETSC_ERROR(err);
    err = PetscFree(tmpVertices);CHECK_PETSC_ERROR(err);
#endif
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while preparing for writing data to HDF5 file "
	<< _filename << " at time " << t << ".\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (const ALE::Exception& err) {
    std::ostringstream msg;
    msg << "Error while preparing for writing data to HDF5 file "
	<< _filename << " at time " << t << ".\n" << err.msg();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Unknown error while preparing for writing data to HDF5 file "
	<< _filename << " at time " << t << ".\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // open

// ----------------------------------------------------------------------
// Close output files.
void
pylith::meshio::DataWriterHDF5Ext<mesh_type,field_type>::close(void)
{ // close
  if (_h5->isOpen()) {
    // :TODO: Update number of time steps.

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
    // :TODO: Must account for possible presence of 'censored depth'
    // and censor the appropriate vertices.
    PetscVec vector = field.vector();
    if (!vector) {
      field.createVector();
      vector = field.vector();
    } // if

#if 0 // TEMPORARY DEBUGGING
    const char* vecname = 0;
    PetscObjectGetName((PetscObject) vector, &vecname);
    std::cout << "NAME field: " << field.label()
	      << ", section: " << field.section()->getName()
	      << ", vec: " << vecname
	      << std::endl;
#endif

    // TODO: Create scatter if necessary
    field.createScatter();
    field.scatterSectionToVector();

    // Create external dataset if necessary
    PetscViewer binaryViewer;
    bool createdExternalDataset = false;
    if (_datasets.find(field.label()) != _datasets.end()) {
      binaryViewer = _datasets[field.label()].viewer;
    } else {
      PetscViewerBinaryCreate(sieveMesh->comm(), &binaryViewer);
      PetscViewerBinaryOpen(sieveMesh->comm(), _datasetFilename(field.label()),
			    FILE_MODE_WRITE, binaryViewer);
      ExternalDataset dataset;
      dataset.numTimeSteps = 0;
      dataset.viewer = binaryViewer;
      _datasets[field.label()] = dataset;

      createdExternalDataset = true;
    } // else

    PetscErrorCode err = VecView(vector, binaryViewer);
    CHECK_PETSC_ERROR(err);
    ++_dataset[field.label()].numTimeSteps;

    if (createdExternalDataset) {
      // Add new external dataset to HDF5 file.
      const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = 
	mesh.sieveMesh();
      assert(!sieveMesh.isNull());
      const std::string labelName = 
	(sieveMesh->hasLabel("censored depth")) ? "censored depth" : "depth";
      const ALE::Obj<typename SieveMesh::numbering_type>& numbering =
	sieveMesh->getFactory()->getNumbering(sieveMesh, labelName, 0);
      assert(!numbering.isNull());

      const ALE::Obj<RealSection>& section = field.section();
      assert(!section.isNull());
      assert(!sieveMesh->getLabelStratum(labelName, 0).isNull());

      const int localFiberDim = 
	(sieveMesh->getLabelStratum(labelName, 0)->size() > 0) ? 
	section->getFiberDimension(*sieveMesh->getLabelStratum(labelName, 0)->begin()) : 0;
      int fiberDim = 0;
      MPI_Allreduce((void *) &localFiberDim, (void *) &fiberDim, 1, 
		    MPI_INT, MPI_MAX, field.mesh().comm());
      assert(fiberDim > 0);

      const int numTimeSteps = DataWriter<mesh_type, field_type>::_numTimeSteps;
      const hsize_t ndims = (numTimeSteps > 0) ? 3 : 2;
      hsize_t* dims = (ndims > 0) ? new int[ndims] : 0;
      if (3 == ndims) {
	dims[0] = 1; // external file only constains 1 time step so far.
	dims[1] = numbering->getGlobalSize()
	dims[2] = fiberDim;
      } else {
	dims[0] = numbering->getGlobalSize()
	dims[1] = fiberDim;
      } // else
      _h5->createDatasetRawExternal("/vertex_fields", field.label(),
				    dims, ndims, H5_NATIVE_DOUBLE);
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
  assert(!_h5);

  try {
    // :TODO: Must account for possible presence of 'censored depth'
    // and censor the appropriate vertices.

    PetscVec vector = field.vector();
    if (vector == PETSC_NULL) {
      field.createVector();
      vector = field.vector();
    }

#if 0 // TEMPORARY DEBUGGING
    const char* vecname = 0;
    PetscObjectGetName((PetscObject) vector, &vecname);
    std::cout << "NAME field: " << field.label()
	      << ", section: " << field.section()->getName()
	      << ", vec: " << vecname
	      << std::endl;
#endif
    // TODO: Create scatter only if necessary
    field.createScatter();
    field.scatterSectionToVector();

    // Create external dataset if necessary
    PetscViewer binaryViewer;
    bool createdExternalDataset = false;
    if (_datasets.find(field.label()) != _datasets.end()) {
      binaryViewer = _datasets[field.label()].viewer;
    } else {
      PetscViewerBinaryCreate(sieveMesh->comm(), &binaryViewer);
      PetscViewerBinaryOpen(sieveMesh->comm(), _datasetFilename(field.label()),
			    FILE_MODE_WRITE, binaryViewer);
      ExternalDataset dataset;
      dataset.numTimeSteps = 0;
      dataset.viewer = binaryViewer;
      _datasets[field.label()] = dataset;

      createdExternalDataset = true;
    } // else

    PetscErrorCode err = VecView(vector, binaryViewer);
    CHECK_PETSC_ERROR(err);
    ++_dataset[field.label()].numTimeSteps;

    if (createdExternalDataset) {
      // Add new external dataset to HDF5 file.
      const ALE::Obj<SieveMesh>& sieveMesh = field.mesh().sieveMesh();
      assert(!sieveMesh.isNull());
      const int cellDepth = (sieveMesh->depth() == -1) ? -1 : 1;
      const int depth = (0 == label) ? cellDepth : labelId;
      const std::string labelName = (0 == label) ?
	((sieveMesh->hasLabel("censored depth")) ?
	 "censored depth" : "depth") : label;
      assert(!sieveMesh->getFactory().isNull());
      const ALE::Obj<typename SieveMesh::numbering_type>& numbering = 
	sieveMesh->getFactory()->getNumbering(sieveMesh, labelName, depth);
      assert(!numbering.isNull());
      assert(!sieveMesh->getLabelStratum(labelName, depth).isNull());
      const ALE::Obj<RealSection>& section = field.section();
      assert(!section.isNull());
      
      const int localFiberDim = 
	(sieveMesh->getLabelStratum(labelName, depth)->size() > 0) ? 
	section->getFiberDimension(*sieveMesh->getLabelStratum(labelName, depth)->begin()) : 0;
      int fiberDim = 0;
      MPI_Allreduce((void *) &localFiberDim, (void *) &fiberDim, 1, 
		    MPI_INT, MPI_MAX, field.mesh().comm());
      assert(fiberDim > 0);

      const int numTimeSteps = DataWriter<mesh_type, field_type>::_numTimeSteps;
      const hsize_t ndims = (numTimeSteps > 0) ? 3 : 2;
      hsize_t* dims = (ndims > 0) ? new int[ndims] : 0;
      if (3 == ndims) {
	dims[0] = 1; // external file only constains 1 time step so far.
	dims[1] = numbering->getGlobalSize()
	dims[2] = fiberDim;
      } else {
	dims[0] = numbering->getGlobalSize()
	dims[1] = fiberDim;
      } // else
      _h5->createDatasetRawExternal("/cell_fields", field.label(),
				    dims, ndims, H5_NATIVE_DOUBLE);
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
  std::ostringstream filename;
  const int indexExt = _filename.find(".h5");
  filename << std::string(_filename, 0, indexExt) << "_" << field << ".h5";

  return std::string(filename.str());
} // _datasetFilename


// End of file 
