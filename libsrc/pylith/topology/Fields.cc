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

#include "Fields.hh" // Implementation of class methods

#include "Field.hh" // USES Field

#include <pylith/utils/error.h> // USES PYLITH_CHECK_ERROR

#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::Fields::Fields(const Mesh& mesh) :
  _mesh(mesh)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::Fields::~Fields(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::Fields::deallocate(void)
{ // deallocate
  const typename map_type::iterator begin = _fields.begin();
  const typename map_type::iterator end = _fields.end();
  for (typename map_type::iterator iter=begin; iter != end; ++iter) {
    delete iter->second; iter->second = 0;
  } // for
} // deallocate

// ----------------------------------------------------------------------
// Check if fields contains a given field.
bool
pylith::topology::Fields::hasField(const char* name) const
{ // hasField
  typename map_type::const_iterator iter = _fields.find(name);
  return iter != _fields.end();
} // hasField

// ----------------------------------------------------------------------
// Add field.
void
pylith::topology::Fields::add(const char* name,
			      const char* label)
{ // add
  if (hasField(name)) {
    std::ostringstream msg;
    msg << "Could not add field '" << name
	<< "' to fields manager, because it already exists.";
    throw std::runtime_error(msg.str());
  } // if
  
  _fields[name] = new Field(_mesh);
  _fields[name]->label(label);
} // add

// ----------------------------------------------------------------------
// Add field.
void 
pylith::topology::Fields::add(const char* name,
			      const char* label,
			      const pylith::topology::FieldBase::DomainEnum domain,
			      const int fiberDim)
{ // add
  if (hasField(name)) {
    std::ostringstream msg;
    msg << "Could not add field '" << name
	<< "' to fields manager, because it already exists.";
    throw std::runtime_error(msg.str());
  } // if
  
  _fields[name] = new Field(_mesh);
  _fields[name]->label(label);
  _fields[name]->newSection(domain, fiberDim);
} // add

// ----------------------------------------------------------------------
// Delete field.
void
pylith::topology::Fields::del(const char* name)
{ // del
  typename map_type::iterator iter = _fields.find(name);
  if (iter == _fields.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name
	<< "' in fields manager to delete.";
    throw std::runtime_error(msg.str());
  } // if
  delete iter->second; iter->second = 0;
  _fields.erase(name);
} // del

// ----------------------------------------------------------------------
// Delete field.
void
pylith::topology::Fields::delField(const char* name)
{ // delField
  del(name);
} // delField

// ----------------------------------------------------------------------
// Get field.
const pylith::topology::Field&
pylith::topology::Fields::get(const char* name) const
{ // get
  typename map_type::const_iterator iter = _fields.find(name);
  if (iter == _fields.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name
	<< "' in fields manager for retrieval.";
    throw std::runtime_error(msg.str());
  } // if
  return *iter->second;
} // get
	   
// ----------------------------------------------------------------------
// Get field.
pylith::topology::Field&
pylith::topology::Fields::get(const char* name)
{ // get
  typename map_type::iterator iter = _fields.find(name);
  if (iter == _fields.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name
	<< "' in fields manager for retrieval.";
    throw std::runtime_error(msg.str());
  } // if
  return *iter->second;
} // get

// ----------------------------------------------------------------------
// Copy layout to other fields.
void
pylith::topology::Fields::copyLayout(const char* name)
{ // copyLayout
  typename map_type::const_iterator src = _fields.find(name);
  if (src == _fields.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name
	<< "' in fields manager for retrieval.";
    throw std::runtime_error(msg.str());
  } // if

  const typename map_type::iterator begin = _fields.begin();
  const typename map_type::iterator end = _fields.end();
  for (typename map_type::iterator iter=begin; iter != end; ++iter)
    if (iter != src)
      iter->second->cloneSection(*src->second);
} // copyLayout

// ----------------------------------------------------------------------
// Get mesh associated with fields.
const pylith::topology::Mesh&
pylith::topology::Fields::mesh(void) const
{ // mesh
  return _mesh;
} // mesh

// ----------------------------------------------------------------------
// Get names of all fields
void
pylith::topology::Fields::fieldNames(int* numNames, 
				     char*** names) const
{ // fieldNames
  assert(numNames);
  assert(names);

  *numNames = _fields.size();
  *names = new char*[_fields.size()];
  assert(*names);
  const typename map_type::const_iterator namesEnd = _fields.end();
  int i = 0;
  for (typename map_type::const_iterator name = _fields.begin(); 
       name != namesEnd;
       ++name) {
    const char len = name->first.length();
    char* newName = 0;
    if (len > 0) {
      newName = new char[len+1];
      strncpy(newName, name->first.c_str(), len+1);
    } else {
      newName = new char[1];
      newName[0] ='\0';
    } // if/else
    (*names)[i++] = newName;
  } // for
} // fieldNames


// End of file 
