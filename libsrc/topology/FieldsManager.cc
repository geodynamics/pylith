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

#include "FieldsManager.hh" // implementation of class methods

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::topology::FieldsManager::FieldsManager(const ALE::Obj<Mesh>& mesh) :
  _mesh(mesh)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::topology::FieldsManager::~FieldsManager(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Add field.
void
pylith::topology::FieldsManager::addReal(const char* name)
{ // addReal
  assert(!_mesh.isNull());

  map_real_type::iterator iter = _real.find(name);
  if (iter != _real.end()) {
    std::ostringstream msg;
    msg << "Could not add field '" << name
	<< "', because it already exists.";
    throw std::runtime_error(msg.str());
  } // if
  
  ALE::Obj<real_section_type> field = 
    new real_section_type(_mesh->comm(), _mesh->debug());
  _real[name] = field;
} // addReal

// ----------------------------------------------------------------------
// Get field.
const ALE::Obj<pylith::real_section_type>&
pylith::topology::FieldsManager::getReal(const char* name)
{ // getReal
  map_real_type::const_iterator iter = _real.find(name);
  if (iter == _real.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if
  return _real[name];
} // getReal


// End of file 
