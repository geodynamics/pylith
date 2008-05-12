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

#include "EventLogger.hh" // Implementation of class methods


#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::utils::EventLogger::EventLogger(void) :
  _className(""),
  _classId(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::utils::EventLogger::~EventLogger(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Setup logging class.
void
pylith::utils::EventLogger::initialize(void)
{ // initialize
  if (_className == "")
    throw std::logic_error("Must set logging class name before "
			   "initializaing EventLogger.");
  
  _events.clear();
  PetscErrorCode err = PetscCookieRegister(_className.c_str(), &_classId);
  if (err) {
    std::ostringstream msg;
    msg << "Could not register logging class '" << _className << "'.";
    throw std::runtime_error(msg.str());
  } // if
  assert(0 != _classId);
} // initialize

// ----------------------------------------------------------------------
// Register event.
int
pylith::utils::EventLogger::registerEvent(const char* name)
{ // registerEvent
  assert(0 != _classId);
  int id = 0;
  PetscErrorCode err = PetscLogEventRegister(name, _classId, &id);
  if (err) {
    std::ostringstream msg;
    msg << "Could not register logging event '" << name
	<< "' for logging class '" << _className << "'.";
    throw std::runtime_error(msg.str());
  } // if  
  _events[name] = id;
  return id;
} // registerEvent

// ----------------------------------------------------------------------
// Get event identifier.
int
pylith::utils::EventLogger::eventId(const char* name)
{ // eventId
  map_event_type::iterator id = _events.find(name);
  if (id == _events.end()) {
    std::ostringstream msg;
    msg << "Could not find logging event '" << name
	<< "' in logging class '" << _className << "'.";
    throw std::runtime_error(msg.str());
  } // if

  return id->second;
} // eventId


// End of file 
