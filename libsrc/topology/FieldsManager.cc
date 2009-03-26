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
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::topology::FieldsManager::FieldsManager(const ALE::Obj<Mesh>& mesh) :
  _mesh(mesh),
  _solutionName("")
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
  
  _real[name] = new real_section_type(_mesh->comm(), _mesh->debug());
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
  return iter->second;
} // getReal

// ----------------------------------------------------------------------
// Remove field.
void
pylith::topology::FieldsManager::delReal(const char* name)
{ // delReal
  map_real_type::const_iterator iter = _real.find(name);
  if (iter == _real.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name << "' to delete.";
    throw std::runtime_error(msg.str());
  } // if
  _real.erase(name);
} // delReal

// ----------------------------------------------------------------------
// Set fiber dimension for field.
void
pylith::topology::FieldsManager::setFiberDimension(const char* name,
						   const int fiberDim,
						   const char* points)
{ // setFiberDimension
  assert(!_mesh.isNull());
  assert(fiberDim >= 0);

  map_real_type::const_iterator iter = _real.find(name);
  if (iter == _real.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name << "' to delete.";
    throw std::runtime_error(msg.str());
  } // if

  assert(!_real[name].isNull());
  if (0 == strcasecmp(points, "vertices")) {
    const ALE::Obj<Mesh::label_sequence>& vertices = _mesh->depthStratum(0);

    if (vertices->size() > 0) {
      _real[name]->setChart(real_section_type::chart_type(*std::min_element(vertices->begin(),vertices->end()),*std::max_element(vertices->begin(),vertices->end())+1));
    } else {
      _real[name]->setChart(real_section_type::chart_type(0,0));
    }
    _real[name]->setFiberDimension(vertices, fiberDim);
  } else if (0 == strcasecmp(points, "cells")) {
    const ALE::Obj<Mesh::label_sequence>& cells = _mesh->heightStratum(0);

    if (cells->size() > 0) {
      _real[name]->setChart(real_section_type::chart_type(*std::min_element(cells->begin(),cells->end()),*std::max_element(cells->begin(),cells->end())+1));
    } else {
      _real[name]->setChart(real_section_type::chart_type(0,0));
    }
    _real[name]->setFiberDimension(cells, fiberDim);
  } else {
    std::ostringstream msg;
    msg << "Could not determine parse '" << points
	<< "' into a known point type when setting fiber dimension to "
	<< fiberDim << " for section '" << name << "'.\n"
	<< "Known point types are 'vertices' and 'cells'.";
    throw std::runtime_error(msg.str());
  } // if/else
} // setFiberDimension

// ----------------------------------------------------------------------
// Allocate field.
void
pylith::topology::FieldsManager::allocate(const char* name)
{ // allocate
  assert(!_mesh.isNull());

  map_real_type::const_iterator iter = _real.find(name);
  if (iter == _real.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name << "' to delete.";
    throw std::runtime_error(msg.str());
  } // if
  
  assert(!_real[name].isNull());
  _mesh->allocate(_real[name]);
} // allocate

// ----------------------------------------------------------------------
// Copy layout of field to all other fields.
void
pylith::topology::FieldsManager::copyLayout(const char* name)
{ // copyLayout
  assert(!_mesh.isNull());

  map_real_type::const_iterator src = _real.find(name);
  if (src == _real.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  assert(!src->second.isNull());  
  const map_real_type::iterator begin = _real.begin();
  const map_real_type::iterator end = _real.end();
  for (map_real_type::iterator iter=begin; iter != end; ++iter)
    if (iter != src) {
      // Make sure fields are same size
      assert(!iter->second.isNull());
      iter->second->setAtlas(src->second->getAtlas());
      iter->second->allocateStorage();
      iter->second->setBC(src->second->getBC());
    } // if
} // copyLayout

// ----------------------------------------------------------------------
// Copy layout of field to managed fields.
void
pylith::topology::FieldsManager::copyLayout(
				    const ALE::Obj<real_section_type>& field)
{ // copyLayout
  assert(!_mesh.isNull());
  assert(!field.isNull());
  
  const map_real_type::iterator begin = _real.begin();
  const map_real_type::iterator end = _real.end();
  for (map_real_type::iterator iter=begin; iter != end; ++iter) {
    // Make sure fields are same size
    assert(!iter->second.isNull());
    iter->second->setAtlas(field->getAtlas());
    iter->second->allocateStorage();
    iter->second->setBC(field->getBC());
  } // for
} // copyLayout

// ----------------------------------------------------------------------
// Set name of solution field.
void
pylith::topology::FieldsManager::solutionField(const char* name)
{ // solutionField
  map_real_type::const_iterator iter = _real.find(name);
  if (iter == _real.end()) {
    std::ostringstream msg;
    msg << "Cannot use unknown field '" << name 
	<< "' when setting name of solution field.";
    throw std::runtime_error(msg.str());
  } // if
  _solutionName = name;
} // solutionField

// ----------------------------------------------------------------------
// Get solution field.
const ALE::Obj<pylith::real_section_type>&
pylith::topology::FieldsManager::getSolution(void)
{ // getSolution
  if (_solutionName == "")
    throw std::runtime_error("Cannot retrieve solution. Name of solution " \
			     "field has not been specified.");
  return getReal(_solutionName.c_str());
} // getSolution

// ----------------------------------------------------------------------
// Create history manager for a subset of the managed fields.
void
pylith::topology::FieldsManager::createHistory(const char** fields,
					       const int size)
{ // createHistory
  if (size > 0 && 0 != fields) {
    _history.resize(size);
    for (int i=0; i < size; ++i) {
      map_real_type::const_iterator iter = _real.find(fields[i]);
      if (iter == _real.end()) {
	std::ostringstream msg;
	msg << "Cannot use unknown field '" << fields[i] 
	    << "' when creating history.";
	throw std::runtime_error(msg.str());
      } // if
      _history[i] = fields[i];
    } // for
  } // if
} // createHistory

// ----------------------------------------------------------------------
// Shift fields in history. Handles to fields are shifted so that the
// most recent values become associated with the second most recent
// item in the history, etc.
void
pylith::topology::FieldsManager::shiftHistory(void)
{ // shiftHistory

  assert(_history.size() > 0);
  const int size = _history.size();
  ALE::Obj<real_section_type> tmp = _real[_history[size-1]];
  tmp.addRef();
  for (int i=size-1; i > 0; --i)
    _real[_history[i]] = _real[_history[i-1]];
  _real[_history[0]] = tmp;

} // shiftHistory

// ----------------------------------------------------------------------
// Get field in history by position.
const ALE::Obj<pylith::real_section_type>&
pylith::topology::FieldsManager::getFieldByHistory(const int index)
{ // getFieldByHistory
  assert(index < _history.size());
  return getReal(_history[index].c_str());
} // getFieldByHistory


// End of file 
