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

#include "pylith/meshio/meshiofwd.hh" // Forward declarations

#include "pylith/utils/arrayfwd.hh" // USES string_vector

#include <string> // HASA std::string

// Forward declarations -------------------------------------------------

// ExodusII ----------------------------------------------------------
class pylith::meshio::ExodusII { // ExodusII
    friend class TestExodusII; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    ExodusII(void);

    /** Open file.
     *
     * @param filename Name of file.
     */
    ExodusII(const char* filename);

    /// Destructor
    ~ExodusII(void);

    /// Deallocate data structures.
    void deallocate(void);

    /** Set filename.
     *
     * @param filename Name of file
     */
    void filename(const char* name);

    /** Get filename.
     *
     * @returns Name of file
     */
    const char* filename(void) const;

    /// Open file.
    void open(void);

    /// Close file.
    void close(void);

    /** Check if file contains dimension.
     *
     * @param[in] name Name of dimension.
     * @param[out] Id of dimension if found, -1 otherwise.
     * @returns True if file contains dimension, false otherwise.
     */
    bool hasDim(const char* name,
                int* id) const;

    /** Check if file contains attribute.
     *
     * @param[in] name Name of attribute.
     * @param[out] Id of attribute if found, -1 otherwise.
     * @returns True if file contains attribute, false otherwise.
     */
    bool hasAtt(const char* name,
                int* id) const;

    /** Check if file contains variable.
     *
     * @param[in] name Name of variable.
     * @param[out] Id of variable if found, -1 otherwise.
     * @returns True if file contains variable, false otherwise.
     */
    bool hasVar(const char* name,
                int* id) const;

    /** Get dimension from file.
     *
     * @param name Name of dimension.
     * @returns Value of dimension.
     */
    int getDim(const char* name) const;

    /** Get values for variable as an array of PylithScalars.
     *
     * @param values Array of values.
     * @param dims Expected dimensions for variable.
     * @param ndims Number of dimension for variable.
     * @param name Name of variable.
     */
    void getVar(PylithScalar* values,
                int* dims,
                int ndims,
                const char* name) const;

    /** Get values for variable as an array of ints.
     *
     * @param values Array of values.
     * @param dims Expected dimensions for variable.
     * @param ndims Number of dimension for variable.
     * @param name Name of variable.
     */
    void getVar(int* values,
                int* dims,
                int ndims,
                const char* name) const;

    /** Get values for variable as an array of strings.
     *
     * @param values Array of values.
     * @param dim Number of values.
     * @param name Name of variable.
     */
    void getVar(string_vector* values,
                int dim,
                const char* name) const;

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    std::string _filename; ///< Name of file
    int _file; ///< ExodusII file

}; // ExodusII

// End of file
