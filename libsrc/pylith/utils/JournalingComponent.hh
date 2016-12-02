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
