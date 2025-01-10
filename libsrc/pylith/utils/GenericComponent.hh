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

#include <string> // HASA std::string

// GenericComponent ----------------------------------------------------------
/** @brief C++ object for managing Pyre journals.
 *
 * Provides object name to facilitate consistent names for Pyre journals for
 * use in error and debugging messages.
 */
class pylith::utils::GenericComponent {
    friend class TestGenericComponent; // unit testing

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /// Constructor
    GenericComponent(void);

    /// Destructor
    ~GenericComponent(void);

    /** Set name of component.
     *
     * @param value Name of component.
     */
    void setName(const char* value);

    /** Get name of component.
     *
     * @returns Name of component.
     */
    const char* getName(void) const;

    // PRIVATE METHODS //////////////////////////////////////////////////////
private:

    std::string _name; ///< Name of component (used in journals).

    // PRIVATE METHODS //////////////////////////////////////////////////////
private:

    GenericComponent(const GenericComponent&); ///< Not implemented
    const GenericComponent& operator=(const GenericComponent&); ///< Not implemented

}; // GenericComponent

// End of file
