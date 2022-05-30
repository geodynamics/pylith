// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

#if !defined(pylith_utils_petscoptions_hh)
#define pylith_utils_petscoptions_hh

// Include directives ---------------------------------------------------
#include "utilsfwd.hh" // forward declarations
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/materials/materialsfwd.hh" // USES Material

#include <map> // HASA std::map
#include <string> // HASA std::string

// ================================================================================================
class pylith::utils::PetscDefaults : public pylith::utils::GenericComponent {
    friend class TestPetscDefaults; // unit testing

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    static const int NONE;
    static const int MONITORS;
    static const int SOLVER;
    static const int PARALLEL;
    static const int TESTING;

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Set default PETSc solver options based on solution field and material.
     *
     * @param[in] solution Solution field for problem.
     * @param[in] material Solution field.
     * @param[in] flags Flags for turning on defaults for PETSc options.
     */
    static
    void set(const pylith::topology::Field& solution,
             const pylith::materials::Material* material,
             const int flags);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    PetscDefaults(void); ///< Not implemented
    PetscDefaults(const PetscDefaults &); ///< Not implemented.
    const PetscDefaults& operator=(const PetscDefaults&); ///< Not implemented

}; // class PetscDefaults

// ================================================================================================
class pylith::utils::PetscOptions : public pylith::utils::GenericComponent {
    friend class TestPetscOptions; // unit testing
    friend class _PetscOptions; // Internal

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    PetscOptions(void);

    /// Destructor
    ~PetscOptions(void);

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

}; // PetscOptions

#endif // pylith_utils_petscoptions_hh

// End of file
