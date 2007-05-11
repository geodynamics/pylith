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

#include <fstream> // USES std::ofstream
#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::SolutionIOVTK::SolutionIOVTK(void) :
  _fout(new std::ofstream)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::SolutionIOVTK::~SolutionIOVTK(void)
{ // destructor
  delete _fout; _fout = 0;
} // destructor  

// ----------------------------------------------------------------------
// Open output files.
void
pylith::meshio::SolutionIOVTK::open(const ALE::Obj<ALE::Mesh>& mesh)
{ // open
  assert(0 != _fout);

  _fout->open(_filename.c_str());
  if (!_fout->is_open() || !_fout->good()) {
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
  _fout->close();
} // close

// ----------------------------------------------------------------------
// Write solution topology to file.
void
pylith::meshio::SolutionIOVTK::writeTopology(const ALE::Obj<ALE::Mesh>& mesh,
			        const spatialdata::geocoords::CoordSys* csMesh)
{ // writeTopology
  try {
    // ADD STUFF HERE

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
    // ADD STUFF HERE
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
