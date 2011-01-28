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
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::open(const mesh_type& mesh,
							   const int numTimeSteps,
							   const char* label,
							   const int labelId)
{ // open
  typedef typename mesh_type::SieveMesh SieveMesh;
  typedef typename mesh_type::SieveMesh::label_sequence label_sequence;
  typedef typename mesh_type::SieveMesh::numbering_type numbering_type;
  typedef typename mesh_type::SieveMesh::sieve_type sieve_type;

  DataWriter<mesh_type, field_type>::open(mesh, numTimeSteps, label, labelId);

  try {
    PetscErrorCode err = 0;
    
    const std::string& filename = _hdf5Filename();

    _timesteps.clear();

    err = PetscViewerCreate(mesh.comm(), &_viewer);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerSetType(_viewer, PETSCVIEWERHDF5);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerFileSetMode(_viewer, FILE_MODE_WRITE);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerFileSetName(_viewer, filename.c_str());
    CHECK_PETSC_ERROR(err);

    const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
    assert(!sieveMesh.isNull());
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
    const char* context = DataWriter<mesh_type, field_type>::_context.c_str();
    topology::Field<mesh_type> coordinates(mesh, coordinatesSection, metadata);
    coordinates.label("vertices");
    ALE::Obj<numbering_type> vNumbering = 
      sieveMesh->hasLabel("censored depth") ?
      sieveMesh->getFactory()->getNumbering(sieveMesh, "censored depth", 0) :
      sieveMesh->getFactory()->getNumbering(sieveMesh, 0);
    assert(!vNumbering.isNull());
    //vNumbering->view("VERTEX NUMBERING");
    coordinates.createScatter(vNumbering, context);
    coordinates.scatterSectionToVector(context);
    PetscVec coordinatesVector = coordinates.vector(context);
    assert(coordinatesVector);
    int blockSize = 1;
    err = VecGetBlockSize(coordinatesVector, &blockSize); // get block size
    CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(coordinatesVector, cs->spaceDim()); // bs for output
    CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PushGroup(_viewer, "/geometry"); CHECK_PETSC_ERROR(err);
    err = VecView(coordinatesVector, _viewer);CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer); CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(coordinatesVector, blockSize); // reset block size
    CHECK_PETSC_ERROR(err);

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
    err = MPI_Allreduce(&numCornersLocal, &numCorners, 1, MPI_INT, MPI_MAX,
		     sieveMesh->comm()); CHECK_PETSC_ERROR(err);

    const int ncells = cNumbering->getLocalSize();
    const int conesSize = ncells*numCorners;
    PetscScalar* tmpVertices = (conesSize > 0) ? new PetscScalar[conesSize] : 0;

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

    PetscVec elemVec;
    err = VecCreateMPIWithArray(sieveMesh->comm(), conesSize, PETSC_DETERMINE,
				tmpVertices, &elemVec); CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) elemVec,
			     "cells");CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PushGroup(_viewer, "/topology"); CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(elemVec, numCorners); CHECK_PETSC_ERROR(err);
    err = VecView(elemVec, _viewer); CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer); CHECK_PETSC_ERROR(err);
    err = VecDestroy(elemVec); CHECK_PETSC_ERROR(err);
    delete[] tmpVertices; tmpVertices = 0;
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while opening HDF5 file "
	<< _filename << ".\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Unknown error while opening HDF5 file "
	<< _filename << ".\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // openTimeStep

// ----------------------------------------------------------------------
// Close output files.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::close(void)
{ // close
  if (_viewer) {
    PetscViewerDestroy(_viewer);
  } // if
  _viewer = 0;
  _timesteps.clear();
} // close

// ----------------------------------------------------------------------
// Write field over vertices to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::writeVertexField(
				            const double t,
					    field_type& field,
					    const mesh_type& mesh)
{ // writeVertexField
  typedef typename mesh_type::SieveMesh::numbering_type numbering_type;

  assert(_viewer);

  try {
    const char* context = DataWriter<mesh_type, field_type>::_context.c_str();

    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
    assert(!sieveMesh.isNull());
    ALE::Obj<numbering_type> vNumbering = 
      sieveMesh->hasLabel("censored depth") ?
      sieveMesh->getFactory()->getNumbering(sieveMesh, "censored depth", 0) :
      sieveMesh->getFactory()->getNumbering(sieveMesh, 0);
    assert(!vNumbering.isNull());
    field.createScatter(vNumbering, context);
    field.scatterSectionToVector(context);
    PetscVec vector = field.vector(context);
    assert(vector);

    const ALE::Obj<typename mesh_type::RealSection>& section = field.section();
    assert(!section.isNull());
    const std::string labelName = 
      (sieveMesh->hasLabel("censored depth")) ? "censored depth" : "depth";
    assert(!sieveMesh->getLabelStratum(labelName, 0).isNull());
    int fiberDimLocal = 
      (sieveMesh->getLabelStratum(labelName, 0)->size() > 0) ? 
      section->getFiberDimension(*sieveMesh->getLabelStratum(labelName, 0)->begin()) : 0;
    int fiberDim = 0;
    MPI_Allreduce(&fiberDimLocal, &fiberDim, 1, MPI_INT, MPI_MAX,
		  field.mesh().comm());
    assert(fiberDim > 0);

    PetscErrorCode err = 0;

    if (_timesteps.find(field.label()) == _timesteps.end())
      _timesteps[field.label()] = 0;
    else
      _timesteps[field.label()] += 1;
    const int istep = _timesteps[field.label()];

    // Set temporary block size that matches fiber dimension for output.
    int blockSize = 0;
    err = VecGetBlockSize(vector, &blockSize); CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(vector, fiberDim); CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PushGroup(_viewer, "/vertex_fields");
    CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5SetTimestep(_viewer, istep); CHECK_PETSC_ERROR(err);
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
  typedef typename mesh_type::SieveMesh::numbering_type numbering_type;

  assert(_viewer);
  
  try {
    const char* context = DataWriter<mesh_type, field_type>::_context.c_str();

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
    field.createScatter(numbering, context);
    field.scatterSectionToVector(context);
    PetscVec vector = field.vector(context);
    assert(vector);

    const ALE::Obj<typename mesh_type::RealSection>& section = field.section();
    assert(!section.isNull());      
    assert(!sieveMesh->getLabelStratum(labelName, depth).isNull());
    int fiberDimLocal = 
      (sieveMesh->getLabelStratum(labelName, depth)->size() > 0) ? 
      section->getFiberDimension(*sieveMesh->getLabelStratum(labelName, depth)->begin()) : 0;
    int fiberDim = 0;
    MPI_Allreduce(&fiberDimLocal, &fiberDim, 1, MPI_INT, MPI_MAX,
		  field.mesh().comm());
    assert(fiberDim > 0);

    PetscErrorCode err = 0;
    
    if (_timesteps.find(field.label()) == _timesteps.end())
      _timesteps[field.label()] = 0;
    else
      _timesteps[field.label()] += 1;
    const int istep = _timesteps[field.label()];

    // Set temporary block size that matches fiber dimension for output.
    int blockSize = 0;
    err = VecGetBlockSize(vector, &blockSize); CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(vector, fiberDim); CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PushGroup(_viewer, "/cell_fields");
    CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5SetTimestep(_viewer, istep); CHECK_PETSC_ERROR(err);
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
