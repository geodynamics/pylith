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
// Copyright (c) 2010-2013 University of California, Davis
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
  DataWriter<mesh_type, field_type>::deallocate();

  PetscErrorCode err = 0;
  err = PetscViewerDestroy(&_viewer); CHECK_PETSC_ERROR(err);
  err = VecDestroy(&_tstamp); CHECK_PETSC_ERROR(err);
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
  typedef typename mesh_type::SieveMesh SieveMesh;
  typedef typename mesh_type::SieveMesh::label_sequence label_sequence;
  typedef typename mesh_type::SieveMesh::numbering_type numbering_type;
  typedef typename mesh_type::SieveMesh::sieve_type sieve_type;

  DataWriter<mesh_type, field_type>::open(mesh, numTimeSteps, label, labelId);

  try {
    PetscErrorCode err = 0;

    deallocate();
    
    const std::string& filename = _hdf5Filename();

    const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
    assert(!sieveMesh.isNull());

    _timesteps.clear();
    _tstampIndex = 0;
    const int commRank = sieveMesh->commRank();
    const int localSize = (!commRank) ? 1 : 0;
    err = VecCreateMPI(mesh.comm(), localSize, 1, &_tstamp);
    CHECK_PETSC_ERROR(err);
    assert(_tstamp);
    err = VecSetBlockSize(_tstamp, 1); CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) _tstamp, "time");

    err = PetscViewerHDF5Open(mesh.comm(), filename.c_str(), FILE_MODE_WRITE,
			      &_viewer);
    CHECK_PETSC_ERROR(err);

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

    coordinates.createScatterWithBC(mesh, vNumbering, context);
    coordinates.scatterSectionToVector(context);
    PetscVec coordinatesVector = coordinates.vector(context);
    assert(coordinatesVector);
    err = PetscViewerHDF5PushGroup(_viewer, "/geometry");CHECK_PETSC_ERROR(err);
    err = VecView(coordinatesVector, _viewer);CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer); CHECK_PETSC_ERROR(err);

    // Account for censored cells
    int cellDepthLocal = (sieveMesh->depth() == -1) ? -1 : 1;
    int cellDepth = 0;
    err = MPI_Allreduce(&cellDepthLocal, &cellDepth, 1, MPI_INT, MPI_MAX, 
			sieveMesh->comm());CHECK_PETSC_ERROR(err);
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
    PylithScalar* tmpVertices = (conesSize > 0) ? new PylithScalar[conesSize] : 0;

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
    err = VecCreateMPIWithArray(sieveMesh->comm(), numCorners, conesSize, PETSC_DETERMINE,
				tmpVertices, &elemVec);CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) elemVec, "cells");CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PushGroup(_viewer, "/topology");CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(elemVec, numCorners);CHECK_PETSC_ERROR(err);
    err = VecView(elemVec, _viewer);CHECK_PETSC_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer);CHECK_PETSC_ERROR(err);
    err = VecDestroy(&elemVec);CHECK_PETSC_ERROR(err);
    delete[] tmpVertices; tmpVertices = 0;

    hid_t h5 = -1;
    err = PetscViewerHDF5GetFileId(_viewer, &h5); CHECK_PETSC_ERROR(err);
    assert(h5 >= 0);
    const int cellDim = mesh.dimension();
    HDF5::writeAttribute(h5, "/topology/cells", "cell_dim", (void*)&cellDim,
			 H5T_NATIVE_INT);
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while opening HDF5 file "
	<< _hdf5Filename() << ".\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Unknown error while opening HDF5 file "
	<< _hdf5Filename() << ".\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // open

// ----------------------------------------------------------------------
// Close output files.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::close(void)
{ // close
  PetscErrorCode err = 0;
  err = PetscViewerDestroy(&_viewer); CHECK_PETSC_ERROR(err);
  err = VecDestroy(&_tstamp); CHECK_PETSC_ERROR(err);
  assert(!_tstamp);

  _timesteps.clear();
  _tstampIndex = 0;

  int rank = 0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    Xdmf metafile;
    const std::string& hdf5filename = _hdf5Filename();
    const int indexExt = hdf5filename.find(".h5");
    std::string xdmfFilename = 
      std::string(hdf5filename, 0, indexExt) + ".xmf";
    metafile.write(xdmfFilename.c_str(), _hdf5Filename().c_str());
  } // if
} // close

// ----------------------------------------------------------------------
// Write field over vertices to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::writeVertexField(
				            const PylithScalar t,
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

#if 0
    std::cout << "WRITE VERTEX FIELD" << std::endl;
    mesh.view("MESH");
    field.view("FIELD");;
    vNumbering->view("NUMBERING");
    std::cout << std::endl;
#endif

    field.createScatterWithBC(mesh, vNumbering, context);
    field.scatterSectionToVector(context);
    PetscVec vector = field.vector(context);
    assert(vector);

    if (_timesteps.find(field.label()) == _timesteps.end())
      _timesteps[field.label()] = 0;
    else
      _timesteps[field.label()] += 1;
    const int istep = _timesteps[field.label()];
    // Add time stamp to "/time" if necessary.
    const int commRank = sieveMesh->commRank();
    if (_tstampIndex == istep)
      _writeTimeStamp(t, commRank);

#if 0 // debugging
    field.view("writeVertexField");
    VecView(vector, PETSC_VIEWER_STDOUT_WORLD);
#endif

    PetscErrorCode err = 0;
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
} // writeVertexField

// ----------------------------------------------------------------------
// Write field over cells to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::writeCellField(
				       const PylithScalar t,
				       field_type& field,
				       const char* label,
				       const int labelId)
{ // writeCellField
  typedef typename mesh_type::SieveMesh::numbering_type numbering_type;

  assert(_viewer);
  
  try {
    const char* context = DataWriter<mesh_type, field_type>::_context.c_str();
    PetscErrorCode err = 0;

    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = 
      field.mesh().sieveMesh();
    assert(!sieveMesh.isNull());
    int cellDepthLocal = (sieveMesh->depth() == -1) ? -1 : 1;
    int cellDepth = 0;
    err = MPI_Allreduce(&cellDepthLocal, &cellDepth, 1, MPI_INT, MPI_MAX, 
			sieveMesh->comm());CHECK_PETSC_ERROR(err);
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

    if (_timesteps.find(field.label()) == _timesteps.end())
      _timesteps[field.label()] = 0;
    else
      _timesteps[field.label()] += 1;
    const int istep = _timesteps[field.label()];
    // Add time stamp to "/time" if necessary.
    const int commRank = sieveMesh->commRank();
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


// ----------------------------------------------------------------------
// Write time stamp to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::_writeTimeStamp(
						    const PylithScalar t,
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
