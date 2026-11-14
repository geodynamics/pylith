// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "petscfwd.hh" // forward declarations

#include <map> // USES std::map
#include <string> // USES std::string

#include "petsc.h"
#include "petsclog.h" // USES PetscLogEventBegin/End() in inline methods

// EventLogger ------------------------------------------------------------------------------------
/** @brief C++ object for managing event logging using PETSc.
 *
 * Each logger object manages the events for a single "logging class".
 */
class pylith::petsc::EventLogger {
    friend class TestEventLogger; // unit testing

    // Scope --------------------------------------------------------------------------------------
  public:
    /** @brief Scope guard that logs an event for the lifetime of the object.
     *
     * The event begins in the constructor and ends in the destructor, so begin/end calls
     * always match, including early returns and exceptions.
     *
     * Usage: `pylith::petsc::EventLogger::Scope logScope(eventId);`
     */
    class Scope {
      public:
        /** Constructor logs event begin.
         *
         * @param[in] event PETSc event identifier from registerEvent().
         */
        explicit Scope(const PetscLogEvent event);

        /// Destructor logs event end.
        ~Scope(void);

      private:
        Scope(const Scope&) = delete;
        Scope& operator=(const Scope&) = delete;

        PetscLogEvent _event; ///< PETSc event identifier.
    }; // Scope

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
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
    PetscLogEvent registerEvent(const char* name);

    /** Register stage.
     *
     * @prerequisite Must call initialize() before registerStage().
     *
     * @param name Name of stage.
     * @returns Stage identifier.
     */
    PetscLogStage registerStage(const char* name);

    /** Get stage identifier.
     *
     * @param name Name of stage.
     * @returns Stage identifier.
     */
    PetscLogStage getStageId(const char* name);

    /** Log stage begin.
     *
     * @param id Stage identifier.
     */
    void stagePush(const PetscLogStage id);

    /// Log stage end.
    void stagePop(void);

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
  private:
    EventLogger(const EventLogger&); ///< Not implemented
    const EventLogger& operator=(const EventLogger&); ///< Not implemented

    // PRIVATE TYPEDEFS ///////////////////////////////////////////////////////////////////////////
  private:
    typedef std::map<std::string, PetscLogEvent> map_event_type;
    typedef std::map<std::string, PetscLogStage> map_stage_type;

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
  private:
    std::string _className; ///< Name of logging class
    PetscClassId _classId; ///< PETSc logging identifier for class
    map_event_type _events; ///< PETSc logging identifiers for events
    map_stage_type _stages; ///< PETSc logging identifiers for stages

}; // EventLogger

#include "EventLogger.icc" // inline methods

// End of file
