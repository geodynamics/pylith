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
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/materials/materialsfwd.hh" // USES Material

#include <map> // HASA std::map
#include <string> // HASA std::string

// ================================================================================================
class pylith::petsc::Defaults : public pylith::utils::GenericComponent {
    friend class TestDefaults; // unit testing

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    static const int NONE;
    static const int MONITORS;
    static const int SOLVER;
    static const int PARALLEL;
    static const int INITIAL_GUESS;
    static const int TESTING;
    static const int COLLECTIVE_IO;
    static const int TS_ADAPTIVE;
    static const int TS_IMPULSE; // :KLUDGE: Flag indicating use of TSAdaptImpulse (full-scale testing)

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Set default PETSc solver options based on solution field and material.
     *
     * @param[in] material Material in problem.
     * @param[in] hasFault True if problem has at least one fault.
     * @param[in] flags Flags for turning on defaults for PETSc options.
     */
    static
    void set(const pylith::materials::Material* material,
             const bool hasFault,
             const int flags);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    Defaults(void) = delete;
    Defaults(const Defaults &) = delete;
    const Defaults& operator=(const Defaults&) = delete;

}; // class Defaults

// ================================================================================================
class pylith::petsc::Options : public pylith::utils::GenericComponent {
    friend class TestOptions; // unit testing
    friend class _Options; // Internal

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    Options(void);

    /// Destructor
    ~Options(void);

    /** Add PETSc option.
     *
     * @param[in] name Option name.
     * @param[in] value Option value.
     */
    void add(const char* name,
             const char* value="");

    /** Remove PETSc option.
     *
     * @param[in] name Option name.
     */
    void remove(const char* name);

    /// Clear PETSc options.
    void clear(void);

    /** Set PETSc options without overriding any previously set options.
     *
     * Options are cleared after being set.
     */
    void set(void);

    /** Set PETSc options, overriding any previously set options with the same name.
     *
     * Options are cleared after being set.
     */
    void override (void);

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    typedef std::map<std::string, std::string> options_t; ///< Map of option name to value.
    options_t _options; ///< Map with PETSc options.

}; // Options

// End of file
