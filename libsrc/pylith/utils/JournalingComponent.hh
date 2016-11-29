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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/utils/JournalingComponent.hh
 *
 * @brief C++ object for managing event logging using PETSc.
 *
 * Each logger object manages the events for a single "logging class".
 */

#if !defined(pylith_utils_journalingcomponent_hh)
#define pylith_utils_journalingcomponent_hh

// Include directives ---------------------------------------------------
#include "utilsfwd.hh" // forward declarations

#include <string> // HASA std::string

// Forward declarations
namespace journal {
    class SeverityDebug;
    class SeverityInfo;
    class SeverityError;
    typedef SeverityDebug debug_t;
    typedef SeverityInfo info_t;
    typedef SeverityError error_t;
} // journal

// JournalingComponent ----------------------------------------------------------
/** @brief C++ object for managing Pyre journals.
 *
 * Provides macros for easy inclusion of Pyre journals when used as a base class.
 */
class pylith::utils::JournalingComponent
{     // JournalingComponent
friend class TestJournalingComponent;     // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

/// Constructor
JournalingComponent(void);

/// Destructor
~JournalingComponent(void);

/** Set name of journal.
 *
 * @param value Name of journal.
 */
void name(const char* value);

/** Get name of journal.
 *
 * @returns Name of journal.
 */
const char* name(void) const;

/// Setup journaling..
void initialize(void);

/** Get debug journal.
 *
 * @returns Debugging journal.
 */
journal::debug_t& debug(void);

/** Get info journal.
 *
 * @returns Debugging journal.
 */
journal::info_t& info(void);

/** Get error journal.
 *
 * @returns Error journal.
 */
journal::error_t& error(void);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected:

journal::debug_t* _debug;     ///< Debugging journal.
journal::info_t* _info;     ///< Info journal.
journal::error_t* _error;     ///< Error journal.

// PRIVATE METHODS //////////////////////////////////////////////////////
private:

std::string _name; ///< Name for journals.

// PRIVATE METHODS //////////////////////////////////////////////////////
private:

JournalingComponent(const JournalingComponent&);     ///< Not implemented
const JournalingComponent& operator=(const JournalingComponent&);     ///< Not implemented

};      // JournalingComponent

#endif // pylith_utils_journalingcomponent_hh


// End of file
