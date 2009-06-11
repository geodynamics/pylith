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

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR

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
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::SolutionFields::deallocate(void)
{ // deallocate
  Fields<Field<Mesh> >::deallocate();
} // deallocate
  
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


// End of file 
