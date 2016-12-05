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
 * @file libsrc/utils/PyreComponent.hh
 *
 * @brief C++ object for managing Pyre component and identifier information.
 */

#if !defined(pylith_utils_pyrecomponent_hh)
#define pylith_utils_pyrecomponent_hh

// Include directives ---------------------------------------------------
#include "utilsfwd.hh" // forward declarations

#include <string> // HASA std::string

// PyreComponent ----------------------------------------------------------
/** @brief C++ object for managing Pyre journals.
 *
 * Provides component name to facilitate consistent names for Pyre journals and
 * component identifier for use in error and debugging messages.
 */
class pylith::utils::PyreComponent
{     // PyreComponent
friend class TestPyreComponent;     // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

/// Constructor
PyreComponent(void);

/// Destructor
~PyreComponent(void);

/** Set name of component.
 *
 * @param value Name of component.
 */
void name(const char* value);

/** Get name of component.
 *
 * @returns Name of component.
 */
const char* name(void) const;

/** Set component identifier (identifies object in component hierarchy).
 *
 * @param value Component identifier.
 */
void identifier(const char* value);

/** Get component identifier (identifies object in component hierarchy).
 *
 * @returns Component identifier.
 */
const char* identifier(void) const;

// PRIVATE METHODS //////////////////////////////////////////////////////
private:

std::string _name; ///< Name of component (used in journals).
std::string _identifier; ///< Identifier for object in component hierarchy.

// PRIVATE METHODS //////////////////////////////////////////////////////
private:

PyreComponent(const PyreComponent&);     ///< Not implemented
const PyreComponent& operator=(const PyreComponent&);     ///< Not implemented

};      // PyreComponent

#endif // pylith_utils_pyrecomponent_hh


// End of file
