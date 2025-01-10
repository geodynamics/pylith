// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/utils/utilsfwd.hh" // forward declarations

#include <string> // USES std::string
#include <map> // USES std::map

#include "petsc.h"
#include "petsclog.h" // USES PetscLogEventBegin/End() in inline methods

// EventLogger ----------------------------------------------------------
/** @brief C++ object for managing event logging using PETSc.
 *
 * Each logger object manages the events for a single "logging class".
 */
class pylith::utils::EventLogger {
    friend class TestEventLogger; // unit testing

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /// Constructor
    EventLogger(void);

    /// Destructor
    ~EventLogger(void);

    /** Set name of logging class.
     *
     * @param name Name of logging class.
     */
    void setClassName(const char* name);

    /** Get name of logging class.
     *
     * @returns Name of logging class.
     */
    const char* getClassName(void) const;

    /// Setup logging class.
    void initialize(void);

    /** Register event.
     *
     * @prerequisite Must call initialize() before registerEvent().
     *
     * @param name Name of event.
     * @returns Event identifier.
     */
    int registerEvent(const char* name);

    /** Get event identifier.
     *
     * @param name Name of event.
     * @returns Event identifier.
     */
    int getEventId(const char* name);

    /** Log event begin.
     *
     * @param id Event identifier.
     */
    void eventBegin(const int id);

    /** Log event end.
     *
     * @param id Event identifier.
     */
    void eventEnd(const int id);

    /** Register stage.
     *
     * @prerequisite Must call initialize() before registerStage().
     *
     * @param name Name of stage.
     * @returns Stage identifier.
     */
    int registerStage(const char* name);

    /** Get stage identifier.
     *
     * @param name Name of stage.
     * @returns Stage identifier.
     */
    int getStageId(const char* name);

    /** Log stage begin.
     *
     * @param id Stage identifier.
     */
    void stagePush(const int id);

    /// Log stage end.
    void stagePop(void);

    // PRIVATE METHODS //////////////////////////////////////////////////////
private:

    EventLogger(const EventLogger&); ///< Not implemented
    const EventLogger& operator=(const EventLogger&); ///< Not implemented

    // PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private:

    typedef std::map<std::string,int> map_event_type;

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    std::string _className; ///< Name of logging class
    int _classId; ///< PETSc logging identifier for class
    map_event_type _events; ///< PETSc logging identifiers for events
    map_event_type _stages; ///< PETSc logging identifiers for stages

}; // EventLogger

#include "EventLogger.icc" // inline methods

// End of file
