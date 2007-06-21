// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/* For now, we write each time step to a different file, so we throw
 * the whole implementation into writeField().
 */

#include <portinfo>

#include "SolutionIOVTK.hh" // implementation of class methods

#include <petscmesh_viewers.hh> // USES VTKViewer
#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::SolutionIOVTK::SolutionIOVTK(void) :
  _viewer(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::SolutionIOVTK::~SolutionIOVTK(void)
{ // destructor
  if (_viewer)
    PetscViewerDestroy(_viewer);
  _viewer = NULL;
} // destructor  

// ----------------------------------------------------------------------
// Open output files.
void
pylith::meshio::SolutionIOVTK::open(const ALE::Obj<ALE::Mesh>& mesh)
{ // open
#if 0
  PetscErrorCode err;

  err = PetscViewerCreate(mesh->comm(), &_viewer);
  err = PetscViewerSetType(_viewer, PETSC_VIEWER_ASCII);
  err = PetscViewerSetFormat(_viewer, PETSC_VIEWER_ASCII_VTK);
  err = PetscViewerFileSetName(_viewer, _filename.c_str());
  if (err) {
    std::ostringstream msg;
    msg << "Could not open VTK file '" << _filename
	<< "' for solution output.\n";
    throw std::runtime_error(msg.str());
  } // if

#endif
} // open

// ----------------------------------------------------------------------
// Close output files.
void
pylith::meshio::SolutionIOVTK::close(void)
{ // close
  if
    (_viewer) PetscViewerDestroy(_viewer);
  _viewer = NULL;
} // close

// ----------------------------------------------------------------------
// Write solution topology to file.
void
pylith::meshio::SolutionIOVTK::writeTopology(const ALE::Obj<ALE::Mesh>& mesh,
			        const spatialdata::geocoords::CoordSys* csMesh)
{ // writeTopology
  #if 0
  try {
    PetscErrorCode err = 0;

    err = VTKViewer::writeHeader(_viewer);
    err = VTKViewer::writeVertices(mesh, _viewer);
    err = VTKViewer::writeElements(mesh, _viewer);
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing topology information to VTK file '"
	<< _filename << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing topology information to VTK file '"
	<< _filename << "'.\n";
    throw std::runtime_error(msg.str());
  } // try/catch
#endif
} // writeTopology

// ----------------------------------------------------------------------
// Write field to file.
void
pylith::meshio::SolutionIOVTK::writeField(
				     const double t,
				     const ALE::Obj<real_section_type>& field,
				     const char* name,
				     const ALE::Obj<ALE::Mesh>& mesh)
{ // writeField

  try {
    PetscErrorCode err;

    std::ostringstream buffer;
    const int indexExt = _filename.find(".vtk");
    buffer << std::string(_filename, 0, indexExt) << "_t" << t << ".vtk";

    err = PetscViewerCreate(mesh->comm(), &_viewer);
    err = PetscViewerSetType(_viewer, PETSC_VIEWER_ASCII);
    err = PetscViewerSetFormat(_viewer, PETSC_VIEWER_ASCII_VTK);
    err = PetscViewerFileSetName(_viewer, buffer.str().c_str());
    if (err) {
      std::ostringstream msg;
      msg << "Could not open VTK file '" << buffer.str()
	  << "' for solution output.\n";
      throw std::runtime_error(msg.str());
    } // if

    err = VTKViewer::writeHeader(_viewer);
    err = VTKViewer::writeVertices(mesh, _viewer);
    err = VTKViewer::writeElements(mesh, _viewer);

    buffer.str("");
    buffer << name << "_t" << t;

    // Now we are enforcing a 3D solution
    //   Perhaps we need to push this argument higher
    err = SectionView_Sieve_Ascii(mesh, field, buffer.str().c_str(), _viewer, 3);
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error while writing field '" << name << "' at time " 
	<< t << " to VTK file '" << _filename << "'.\n" << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) { 
    std::ostringstream msg;
    msg << "Error while writing field '" << name << "' at time " 
	<< t << " to VTK file '" << _filename << "'.\n";
    throw std::runtime_error(msg.str());
  } // try/catch
} // writeField


// End of file 
