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

#include "ParameterManager.hh" // implementation of class methods

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ParameterManager::ParameterManager(
					  const ALE::Obj<Mesh>& mesh) :
  _mesh(mesh)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ParameterManager::~ParameterManager(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Add parameter.
void
pylith::feassemble::ParameterManager::addReal(const char* name)
{ // addReal
  map_real_type::iterator iter = _real.find(name);
  if (iter != _real.end()) {
    std::ostringstream msg;
    msg << "Could not add parameter '" << name
	<< "'. Parameter already exists.";
    throw std::runtime_error(msg.str());
  } // if
  
  ALE::Obj<real_section_type> parameter = 
    new real_section_type(_mesh->comm(), _mesh->debug());
  _real[name] = parameter;
} // addReal

// ----------------------------------------------------------------------
// Get parameter.
const ALE::Obj<pylith::feassemble::ParameterManager::Mesh::real_section_type>&
pylith::feassemble::ParameterManager::getReal(const char* name)
{ // getReal
  map_real_type::const_iterator iter = _real.find(name);
  if (iter == _real.end()) {
    std::ostringstream msg;
    msg << "Could not find parameter '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if
  return _real[name];
} // getReal


// End of file 
