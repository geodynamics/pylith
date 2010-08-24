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

#include <petscmesh_viewers.hh> // USES HDF5Viewer

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
#if 0
  // MATT - This stuff needs to be updated for HDF5.

  try {
    PetscErrorCode err = 0;
    
    const std::string& filename = _hdf5Filename(t);

    err = PetscViewerCreate(mesh.comm(), &_viewer);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerSetType(_viewer, PETSCVIEWERASCII);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerSetFormat(_viewer, PETSC_VIEWER_ASCII_HDF5);
    CHECK_PETSC_ERROR(err);
    err = PetscViewerFileSetName(_viewer, filename.c_str());
    CHECK_PETSC_ERROR(err);

    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
    
    err = HDF5Viewer::writeHeader(sieveMesh, _viewer);
    CHECK_PETSC_ERROR(err);
    //std::cout << "Wrote header for " << filename << std::endl;
    err = HDF5Viewer::writeVertices(sieveMesh, _viewer);
    CHECK_PETSC_ERROR(err);
    //std::cout << "Wrote vertices for " << filename << std::endl;
    if (0 == label) {
      err = HDF5Viewer::writeElements(sieveMesh, _viewer);
      CHECK_PETSC_ERROR(err);
    } else {
      const std::string labelName = 
	(sieveMesh->hasLabel("censored depth")) ? "censored depth" : "depth";
      err = HDF5Viewer::writeElements(sieveMesh, label, labelId, labelName, 0, _viewer);      
      CHECK_PETSC_ERROR(err);
    } // if
    //std::cout << "Wrote elements for " << filename << std::endl;

    _wroteVertexHeader = false;
    _wroteCellHeader = false;
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

#endif
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
					    const field_type& field,
					    const mesh_type& mesh)
{ // writeVertexField
#if 0
  // MATT - This stuff needs to be update for HDF5.

  typedef typename mesh_type::SieveMesh SieveMesh;
  typedef typename field_type::Mesh::RealSection RealSection;

  try {
    int rank = 0;
    MPI_Comm_rank(field.mesh().comm(), &rank);

    const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
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
    const int enforceDim =
      (field.vectorFieldType() != topology::FieldBase::VECTOR) ? fiberDim : 3;

    PetscErrorCode err = 0;
    if (!_wroteVertexHeader) {
      err = PetscViewerASCIIPrintf(_viewer, "POINT_DATA %d\n", 
				   numbering->getGlobalSize());
      CHECK_PETSC_ERROR(err);
      _wroteVertexHeader = true;
    } // if

    err = HDF5Viewer::writeField(section, field.label(), fiberDim, numbering,
				_viewer, enforceDim, _precision);
    CHECK_PETSC_ERROR(err);
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

#endif
} // writeVertexField

// ----------------------------------------------------------------------
// Write field over cells to file.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::writeCellField(
				       const double t,
				       const field_type& field,
				       const char* label,
				       const int labelId)
{ // writeCellField
#if 0
  // MATT - This stuff needs to be update to HDF5.

  typedef typename field_type::Mesh::SieveMesh SieveMesh;
  typedef typename field_type::Mesh::RealSection RealSection;

  try {
    int rank = 0;
    MPI_Comm_rank(field.mesh().comm(), &rank);

    // Correctly handle boundary and fault meshes
    //   Cannot just use mesh->depth() because boundaries report the wrong thing
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
    const int enforceDim =
      (field.vectorFieldType() != topology::FieldBase::VECTOR) ? fiberDim : 3;

    PetscErrorCode err = 0;
    if (!_wroteCellHeader) {
      err = PetscViewerASCIIPrintf(_viewer, "CELL_DATA %d\n", 
				   numbering->getGlobalSize());
      CHECK_PETSC_ERROR(err);
      _wroteCellHeader = true;
    } // if

    HDF5Viewer::writeField(section, field.label(), fiberDim, numbering,
			  _viewer, enforceDim, _precision);
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

#endif
} // writeCellField

// ----------------------------------------------------------------------
// Generate filename for HDF5 file.
template<typename mesh_type, typename field_type>
std::string
pylith::meshio::DataWriterHDF5<mesh_type,field_type>::_hdf5Filename(const double t) const
{ // _hdf5Filename
  std::ostringstream filename;
  const int indexExt = _filename.find(".h5");
  const int numTimeSteps = DataWriter<mesh_type, field_type>::_numTimeSteps;
  if (0 == numTimeSteps)
    filename << std::string(_filename, 0, indexExt) << "_info.hdf5";

  return std::string(filename.str());
} // _hdf5Filename


// End of file 
