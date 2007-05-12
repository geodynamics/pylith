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

#include <portinfo>

#include "SolutionIOVTK.hh" // implementation of class methods

#include <src/dm/mesh/meshvtk.h> // USES VTKViewer
#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

template<typename Bundle, typename Section>
PetscErrorCode SectionView_Sieve_Ascii(const Obj<Bundle>& bundle, const Obj<Section>& s, const char name[], PetscViewer viewer);

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::SolutionIOVTK::SolutionIOVTK(void) :
  _viewer(NULL)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::SolutionIOVTK::~SolutionIOVTK(void)
{ // destructor
  if (_viewer) PetscViewerDestroy(_viewer); _viewer = NULL;
} // destructor  

// ----------------------------------------------------------------------
// Open output files.
void
pylith::meshio::SolutionIOVTK::open(const ALE::Obj<ALE::Mesh>& mesh)
{ // open
  PetscErrorCode ierr;

  ierr = PetscViewerCreate(mesh->comm(), &_viewer);
  ierr = PetscViewerSetType(_viewer, PETSC_VIEWER_ASCII);
  ierr = PetscViewerSetFormat(_viewer, PETSC_VIEWER_ASCII_VTK);
  ierr = PetscViewerFileSetName(_viewer, _filename.c_str());
  if (ierr) {
    std::ostringstream msg;
    msg << "Could not open VTK file '" << _filename
	<< "' for solution output.\n";
    throw std::runtime_error(msg.str());
  } // if

  // Write header
  // ADD STUFF HERE
} // open

// ----------------------------------------------------------------------
// Close output files.
void
pylith::meshio::SolutionIOVTK::close(void)
{ // close
  if (_viewer) PetscViewerDestroy(_viewer); _viewer = NULL;
} // close

// ----------------------------------------------------------------------
// Write solution topology to file.
void
pylith::meshio::SolutionIOVTK::writeTopology(const ALE::Obj<ALE::Mesh>& mesh,
			        const spatialdata::geocoords::CoordSys* csMesh)
{ // writeTopology
  try {
    PetscErrorCode ierr;

    ierr = VTKViewer::writeHeader(_viewer);
    ierr = VTKViewer::writeVertices(mesh, _viewer);
    ierr = VTKViewer::writeElements(mesh, _viewer);
    // Use spatialdata::geocoords::Converter::convert() to convert
    // coordinates of vertices from csMesh to _cs (postpone and wait
    // for more general implementation of SolutionIO?).
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
    PetscErrorCode ierr;

    // Ignore time for now
    ierr = SectionView_Sieve_Ascii(mesh, field, name, _viewer);
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
