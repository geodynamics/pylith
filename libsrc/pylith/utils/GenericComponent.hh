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

/**
 * @file libsrc/utils/GenericComponent.hh
 *
 * @brief C++ object for managing generic component and identifier information.
 */

#if !defined(pylith_utils_genericcomponent_hh)
#define pylith_utils_genericcomponent_hh

// Include directives ---------------------------------------------------
#include "utilsfwd.hh" // forward declarations

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

#endif // pylith_utils_genericcomponent_hh

// End of file
