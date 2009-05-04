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

#include "Metadata.hh" // implementation of class methods

#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor.
pylith::materials::Metadata::Metadata(const ParamDescription* props,
				      const int numProps,
				      const char* dbProps[],
				      const int numDBProps,
				      const ParamDescription* vars,
				      const int numVars,
				      const char* dbVars[],
				      const int numDBVars) :
  _numDBProperties(numDBProps),
  _dbProperties(dbProps),
  _numDBStateVars(numDBVars),
  _dbStateVars(dbVars)
{ // constructor
  ParameterInfo info;

  // Set physical property information.
  _properties.clear();
  _propertyNames.resize(numProps);
  for (int i=0; i < numProps; ++i) {
    info.fiberDim = props[i].fiberDim;
    info.fieldType = props[i].fieldType;
    _properties[props[i].name] = info;
    _propertyNames[i] = props[i].name;
  } // for

  // Set state variable information.
  _stateVars.clear();
  _stateVarNames.resize(numVars);
  for (int i=0; i < numVars; ++i) {
    info.fiberDim = vars[i].fiberDim;
    info.fieldType = vars[i].fieldType;
    _stateVars[vars[i].name] = info;
    _stateVarNames[i] = vars[i].name;
  } // for
} // constructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::materials::Metadata::Metadata(const Metadata& m) :
  _properties(m._properties),
  _stateVars(m._stateVars),
  _propertyNames(m._propertyNames),
  _stateVarNames(m._stateVarNames),
  _dbProperties(m._dbProperties),
  _dbStateVars(m._dbStateVars),
  _numDBProperties(m._numDBProperties),
  _numDBStateVars(m._numDBStateVars)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Default destructor
pylith::materials::Metadata::~Metadata(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get fiber dimension of value.
int
pylith::materials::Metadata::fiberDim(const char* name,
				      const ValueEnum valueType) const
{ // fiberDim
  int fiberDim = 0;

  switch(valueType)
    { // switch

    case PROPERTY : {
      ParameterMap::const_iterator iter = _properties.find(name);
      if (iter == _properties.end()) {
	std::ostringstream msg;
	msg << "Could not find property '" << name
	    << "' in list of properties to get fiber dimension.";
	throw std::runtime_error(msg.str());
      } // if
      fiberDim = iter->second.fiberDim;
      break;
    } // PROPERTY

    case STATEVAR : {
      ParameterMap::const_iterator iter = _stateVars.find(name);
      if (iter == _stateVars.end()) {
	std::ostringstream msg;
	msg << "Could not find state variable '" << name
	    << "' in list of state variables to get fiber dimension.";
	throw std::runtime_error(msg.str());
      } // if
      fiberDim = iter->second.fiberDim;
      break;
    } // STATEVAR

    default :
      assert(0);
    } // switch

  return fiberDim;
} // fiberDim

// ----------------------------------------------------------------------
// Get type of vector field associated with value.
pylith::topology::FieldBase::VectorFieldEnum
pylith::materials::Metadata::fieldType(const char* name,
				       const ValueEnum valueType) const
{ // fieldType
  topology::FieldBase::VectorFieldEnum fieldType = topology::FieldBase::OTHER;

  switch(valueType)
    { // switch

    case PROPERTY : {
      ParameterMap::const_iterator iter = _properties.find(name);
      if (iter == _properties.end()) {
	std::ostringstream msg;
	msg << "Could not find property '" << name
	    << "' in list of properties to get vector field type.";
	throw std::runtime_error(msg.str());
      } // if
      fieldType = iter->second.fieldType;
      break;
    } // PROPERTY

    case STATEVAR : {
      ParameterMap::const_iterator iter = _stateVars.find(name);
      if (iter == _stateVars.end()) {
	std::ostringstream msg;
	msg << "Could not find state variable '" << name
	    << "' in list of state variables to get vector field type.";
	throw std::runtime_error(msg.str());
      } // if
      fieldType = iter->second.fieldType;
      break;
    } // STATEVAR
      
    default :
      assert(0);
    } // switch
  
  return fieldType;
} // fieldType
   

// End of file
