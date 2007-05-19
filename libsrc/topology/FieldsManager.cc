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
  const real_section_type::chart_type& srcChart = src->second->getChart();
  const real_section_type::chart_type::iterator srcBegin = srcChart.begin();
  const real_section_type::chart_type::iterator srcEnd = srcChart.end();
  
  const map_real_type::iterator begin = _real.begin();
  const map_real_type::iterator end = _real.end();
  for (map_real_type::iterator iter=begin; iter != end; ++iter)
    if (iter != src) {
      // Make sure fields are same size
      assert(!iter->second.isNull());
      for (real_section_type::chart_type::iterator p_iter=srcBegin;
	   p_iter != srcEnd;
	   ++p_iter) {
	iter->second->setFiberDimension(*p_iter, 
					src->second->getFiberDimension(*p_iter));
	iter->second->setConstraintDimension(*p_iter, 
					     src->second->getConstraintDimension(*p_iter));
      } // for
      _mesh->allocate(iter->second);
      for (real_section_type::chart_type::iterator p_iter=srcBegin;
	   p_iter != srcEnd;
	   ++p_iter)
	iter->second->setConstraintDof(*p_iter, 
				       src->second->getConstraintDof(*p_iter));
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
  
  const real_section_type::chart_type& srcChart = field->getChart();
  const real_section_type::chart_type::iterator srcBegin = srcChart.begin();
  const real_section_type::chart_type::iterator srcEnd = srcChart.end();
  
  const map_real_type::iterator begin = _real.begin();
  const map_real_type::iterator end = _real.end();
  for (map_real_type::iterator iter=begin; iter != end; ++iter) {
    // Make sure fields are same size
    assert(!iter->second.isNull());
    for (real_section_type::chart_type::iterator p_iter=srcBegin;
	 p_iter != srcEnd;
	 ++p_iter) {
      iter->second->setFiberDimension(*p_iter, 
				      field->getFiberDimension(*p_iter));
      iter->second->setConstraintDimension(*p_iter, 
					   field->getConstraintDimension(*p_iter));
    } // for
    _mesh->allocate(iter->second);
    for (real_section_type::chart_type::iterator p_iter=srcBegin;
	 p_iter != srcEnd;
	 ++p_iter)
      iter->second->setConstraintDof(*p_iter, 
				     field->getConstraintDof(*p_iter));
  } // for
} // copyLayout


// End of file 
