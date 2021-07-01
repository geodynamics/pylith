// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "EventLogger.hh" // Implementation of class methods

#include "error.hh" // USES PYLITH_METHOD_BEGIN/END

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::utils::EventLogger::EventLogger(void) :
    _className(""),
    _classId(0) {}


// ----------------------------------------------------------------------
// Destructor
pylith::utils::EventLogger::~EventLogger(void) {}


// ----------------------------------------------------------------------
// Setup logging class.
void
pylith::utils::EventLogger::initialize(void) {
    PYLITH_METHOD_BEGIN;

    if (_className == "") {
        throw std::logic_error("Must set logging class name before initializing EventLogger.");
    }

    _events.clear();
    PetscErrorCode err = PetscClassIdRegister(_className.c_str(), &_classId);
    if (err) {
        std::ostringstream msg;
        msg << "Could not register logging class '" << _className << "'.";
        throw std::runtime_error(msg.str());
    } // if
    assert(_classId);

    PYLITH_METHOD_END;
} // initialize


// ----------------------------------------------------------------------
// Register event.
int
pylith::utils::EventLogger::registerEvent(const char* name) {
    PYLITH_METHOD_BEGIN;

    assert(_classId);
    int id = 0;
    PetscErrorCode err = PetscLogEventRegister(name, _classId, &id);
    if (err) {
        std::ostringstream msg;
        msg << "Could not register logging event '" << name << "' for logging class '" << _className << "'.";
        throw std::runtime_error(msg.str());
    } // if
    _events[name] = id;
    PYLITH_METHOD_RETURN(id);
} // registerEvent


// ----------------------------------------------------------------------
// Get event identifier.
int
pylith::utils::EventLogger::getEventId(const char* name) {
    PYLITH_METHOD_BEGIN;

    map_event_type::iterator iter = _events.find(name);
    if (iter == _events.end()) {
        std::ostringstream msg;
        msg << "Could not find logging event '" << name << "' in logging class '" << _className << "'.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_RETURN(iter->second);
} // getEventId


// ----------------------------------------------------------------------
// Register stage.
int
pylith::utils::EventLogger::registerStage(const char* name) {
    PYLITH_METHOD_BEGIN;

    assert(_classId);
    int id = 0;
    PetscErrorCode err = PetscLogStageRegister(name, &id);
    if (err) {
        std::ostringstream msg;
        msg << "Could not register logging stage '" << name << "'.";
        throw std::runtime_error(msg.str());
    } // if
    _stages[name] = id;

    PYLITH_METHOD_RETURN(id);
} // registerStage


// ----------------------------------------------------------------------
// Get stage identifier.
int
pylith::utils::EventLogger::getStageId(const char* name) {
    PYLITH_METHOD_BEGIN;

    map_event_type::iterator iter = _stages.find(name);
    if (iter == _stages.end()) {
        registerStage(name);
        iter = _stages.find(name);
        assert(iter != _stages.end());
    } // if

    PYLITH_METHOD_RETURN(iter->second);
} // getStageId


// End of file
