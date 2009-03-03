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

#include "SolutionFields.hh" // implementation of class methods

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::SolutionFields::SolutionFields(const Mesh& mesh) :
  Fields<Field<Mesh> >(mesh),
  _solutionName("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::SolutionFields::~SolutionFields(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set name of solution field.
void
pylith::topology::SolutionFields::solutionName(const char* name)
{ // solutionName
  map_type::const_iterator iter = _fields.find(name);
  if (iter == _fields.end()) {
    std::ostringstream msg;
    msg << "Cannot use unknown field '" << name 
	<< "' when setting name of solution field.";
    throw std::runtime_error(msg.str());
  } // if
  _solutionName = name;
} // solutionName

// ----------------------------------------------------------------------
// Get solution field.
const pylith::topology::Field<pylith::topology::Mesh>&
pylith::topology::SolutionFields::solution(void) const
{ // solution
  if (_solutionName == "")
    throw std::runtime_error("Cannot retrieve solution. Name of solution " \
			     "field has not been specified.");
  return get(_solutionName.c_str());
} // solution

// ----------------------------------------------------------------------
// Get solution field.
pylith::topology::Field<pylith::topology::Mesh>&
pylith::topology::SolutionFields::solution(void)
{ // solution
  if (_solutionName == "")
    throw std::runtime_error("Cannot retrieve solution. Name of solution " \
			     "field has not been specified.");
  return get(_solutionName.c_str());
} // solution

// ----------------------------------------------------------------------
// Create history manager for a subset of the managed fields.
void
pylith::topology::SolutionFields::createHistory(const char** fields,
						const int size)
{ // createHistory
  if (size > 0 && 0 != fields) {
    _history.resize(size);
    for (int i=0; i < size; ++i) {
      map_type::const_iterator iter = _fields.find(fields[i]);
      if (iter == _fields.end()) {
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
// Shift fields in history.
void
pylith::topology::SolutionFields::shiftHistory(void)
{ // shiftHistory
  assert(_history.size() > 0);
  const int size = _history.size();
  Field<Mesh>* tmp = _fields[_history[size-1]];
  for (int i=size-1; i > 0; --i)
    _fields[_history[i]] = _fields[_history[i-1]];
  _fields[_history[0]] = tmp;
} // shiftHistory


// End of file 
