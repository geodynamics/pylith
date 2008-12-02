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
#include <cassert> // USES assert()

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
  map_event_type::iterator iter = _events.find(name);
  if (iter == _events.end()) {
    std::ostringstream msg;
    msg << "Could not find logging event '" << name
	<< "' in logging class '" << _className << "'.";
    throw std::runtime_error(msg.str());
  } // if

  return iter->second;
} // eventId

// ----------------------------------------------------------------------
// Register stage.
int
pylith::utils::EventLogger::registerStage(const char* name)
{ // registerStage
  assert(0 != _classId);
  int id = 0;
  PetscErrorCode err = PetscLogStageRegister(name, &id);
  if (err) {
    std::ostringstream msg;
    msg << "Could not register logging stage '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if  
  _stages[name] = id;

  return id;
} // registerStage

// ----------------------------------------------------------------------
// Get stage identifier.
int
pylith::utils::EventLogger::stageId(const char* name)
{ // stageId
  map_event_type::iterator iter = _stages.find(name);
  if (iter == _stages.end()) {
    std::ostringstream msg;
    msg << "Could not find logging stage '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if

  return iter->second;
} // stagesId


// End of file 
