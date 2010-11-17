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

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::DataWriterHDF5(void) :
  _filename("output.h5"),
  _viewer(0)
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
  if (0 != _viewer)
    PetscViewerDestroy(_viewer);
  _viewer = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor.
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::DataWriterHDF5(const DataWriterHDF5<mesh_type, field_type>& w) :
  DataWriter<mesh_type, field_type>(w),
  _filename(w._filename),
  _viewer(0)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Prepare file for data at a new time step.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::openTimeStep(const double t,
						       const mesh_type& mesh,
						       const char* label,
						       const int labelId)
{ // openTimeStep
  try {
    PetscErrorCode err = 0;
    
    const std::string& filename = _hdf5Filename();

    err = PetscViewerCreate(mesh.comm(), &_viewer);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerSetType(_viewer, PETSCVIEWERHDF5);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerFileSetMode(_viewer, FILE_MODE_WRITE);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerFileSetName(_viewer, filename.c_str());
    CHECK_PETSC_ERROR(err);

    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
    const ALE::Obj<typename mesh_type::RealSection>& coordinatesSection = 
      sieveMesh->getRealSection("coordinates");
    const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
    assert(cs);
    topology::FieldBase::Metadata metadata;
    // :KLUDGE: We would like to use field_type for the coordinates
    // field. However, the mesh coordinates are Field<mesh_type> and
    // field_type can be Field<Mesh> (e.g., displacement field over a
    // SubMesh).
    const char* context = DataWriter<mesh_type, field_type>::_context.c_str();
    topology::Field<mesh_type> coordinates(mesh, coordinatesSection, metadata);
    coordinates.label("vertices");
    coordinates.createVector(context);
    coordinates.createScatter(context);
    coordinates.scatterSectionToVector(context);
    const PetscVec coordinatesVector = coordinates.vector(context);
    int blockSize = 1;
    err = VecGetBlockSize(coordinatesVector, &blockSize);
    CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(coordinatesVector, cs->spaceDim());
    CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PushGroup(_viewer, "/geometry"); CHECK_PETSC_ERROR(err);
    err = VecView(coordinatesVector, _viewer);CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer); CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(coordinatesVector, blockSize); // reset
    CHECK_PETSC_ERROR(err);

    Vec          elemVec;
    PetscScalar *tmpVertices;
    PetscBool    columnMajor = PETSC_FALSE;

    // :TODO: Update this to use sizes from numbering to account for
    // censored vertices.
    const ALE::Obj<typename mesh_type::SieveMesh::label_sequence>& cells =
      sieveMesh->heightStratum(0);
    int numCornersLocal = 0;
    if (cells->size() > 0)
      numCornersLocal = sieveMesh->getNumCellCorners(*cells->begin());
    int numCorners = numCornersLocal;
    err = MPI_Reduce(&numCornersLocal, &numCorners, 1, MPI_INT, MPI_MAX, 0, 
		     sieveMesh->comm()); CHECK_PETSC_ERROR(err);

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
    err = PetscViewerHDF5PushGroup(_viewer, "/topology"); CHECK_PETSC_ERROR(err);
#if 0 // :TODO: ONLY WORKS IF CELLS ARE CENSORED
    err = VecSetBlockSize(elemVec, numCorners); CHECK_PETSC_ERROR(err);
#endif
    err = VecView(elemVec, _viewer);CHECK_PETSC_ERROR(err);
    err = VecDestroy(elemVec);CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer); CHECK_PETSC_ERROR(err);
    err = PetscFree(tmpVertices);CHECK_PETSC_ERROR(err);
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while preparing for writing data to HDF5 file "
	<< _filename << " at time " << t << ".\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Unknown error while preparing for writing data to HDF5 file "
	<< _filename << " at time " << t << ".\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // openTimeStep

// ----------------------------------------------------------------------
/// Cleanup after writing data for a time step.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::closeTimeStep(void)
{ // closeTimeStep
  PetscViewerDestroy(_viewer); _viewer = 0;
} // closeTimeStep

// ----------------------------------------------------------------------
// Write field over vertices to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::writeVertexField(
				            const double t,
					    field_type& field,
					    const mesh_type& mesh)
{ // writeVertexField
  try {
    // We will try the simplest thing, using the embedded vector. If this is not
    // general enough, due to ordering, etc., we can construct an auxiliary vector.

    const char* context = DataWriter<mesh_type, field_type>::_context.c_str();
    field.createVector(context);
    PetscVec vector = field.vector(context);
    assert(vector);

#if 0 // TEMPORARY DEBUGGING
    const char* vecname = 0;
    PetscObjectGetName((PetscObject) vector, &vecname);
    std::cout << "NAME field: " << field.label()
	      << ", section: " << field.section()->getName()
	      << ", vec: " << vecname
	      << std::endl;
#endif
    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
    assert(!sieveMesh.isNull());

    if (sieveMesh->hasLabel("censored depth")) {
      // Remove Lagrange vertices
      const Obj<typename Mesh::numbering_type> vNumbering = sieveMesh->getFactory()->getNumbering(sieveMesh, "censored depth", 0);

      field.createScatter(vNumbering, context);
    } else {
      field.createScatter(context);
    } // if/else
    field.scatterSectionToVector(context);

    const ALE::Obj<typename mesh_type::RealSection>& section = field.section();
    assert(!section.isNull());
    const std::string labelName = 
      (sieveMesh->hasLabel("censored depth")) ? "censored depth" : "depth";
    assert(!sieveMesh->getLabelStratum(labelName, 0).isNull());
    const int fiberDimLocal = 
      (sieveMesh->getLabelStratum(labelName, 0)->size() > 0) ? 
      section->getFiberDimension(*sieveMesh->getLabelStratum(labelName, 0)->begin()) : 0;
    int fiberDim = 0;
    MPI_Allreduce((void *) &fiberDimLocal, (void *) &fiberDim, 1, 
		  MPI_INT, MPI_MAX, field.mesh().comm());
    assert(fiberDim > 0);

    PetscErrorCode err = 0;
    err = PetscViewerHDF5PushGroup(_viewer, "/vertex_fields"); CHECK_PETSC_ERROR(err);
    int blockSize = 0;
    err = VecGetBlockSize(vector, &blockSize); CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(vector, fiberDim); CHECK_PETSC_ERROR(err);
    err = VecView(vector, _viewer); CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer); CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(vector, blockSize); CHECK_PETSC_ERROR(err);
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to HDF5 file '" << _filename << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to HDF5 file '" << _filename << "'.\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // writeVertexField

// ----------------------------------------------------------------------
// Write field over cells to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::writeCellField(
				       const double t,
				       field_type& field,
				       const char* label,
				       const int labelId)
{ // writeCellField
  try {
    // We will try the simplest thing, using the embedded vector. If this is not
    // general enough, due to ordering, etc., we can construct an auxiliary vector.

    const char* context = DataWriter<mesh_type, field_type>::_context.c_str();
      field.createVector(context);
    PetscVec vector = field.vector(context);
    assert(vector);

#if 0 // TEMPORARY DEBUGGING
    const char* vecname = 0;
    PetscObjectGetName((PetscObject) vector, &vecname);
    std::cout << "NAME field: " << field.label()
	      << ", section: " << field.section()->getName()
	      << ", vec: " << vecname
	      << std::endl;
#endif
    field.createScatter(context);
    field.scatterSectionToVector(context);

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

    PetscErrorCode err = 0;
    err = PetscViewerHDF5PushGroup(_viewer, "/cell_fields"); CHECK_PETSC_ERROR(err);
    
    int blockSize = 0;
    err = VecGetBlockSize(vector, &blockSize); CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(vector, fiberDim); CHECK_PETSC_ERROR(err);
    err = VecView(vector, _viewer); CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer); CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(vector, blockSize); CHECK_PETSC_ERROR(err);
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to HDF5 file '" << _filename << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << field.label() << "' at time " 
	<< t << " to HDF5 file '" << _filename << "'.\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // writeCellField

// ----------------------------------------------------------------------
// Generate filename for HDF5 file.
template<typename mesh_type, typename field_type>
std::string
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::_hdf5Filename(void) const
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


// End of file 
