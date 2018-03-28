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
 * @file libsrc/meshio/OutputSoln.hh
 *
 * @brief C++ object for managing output of the solution over the domain.
 */

#if !defined(pylith_meshio_OutputSoln_hh)
#define pylith_meshio_OutputSoln_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "OutputManager.hh" // ISA OutputManager

#include "pylith/utils/array.hh" // HASA string_vector

// OutputSoln -----------------------------------------------------
/** @brief C++ object for managing output of the solution over the domain.
 */
class pylith::meshio::OutputSoln : public OutputManager {
    friend class TestOutputSoln;   // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    OutputSoln(void);

    /// Destructor
    ~OutputSoln(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set names of solution fields to output.
     *
     * @param[in] names Array of names of fields to output.
     * @param[in] numNames Length of array.
     */
    void vertexDataFields(const char* names[],
                          const int numNames);

    /** Verify configuration.
     *
     * @param mesh PETSc mesh
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Write solution at time step.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     */
    void writeTimeStep(const PylithReal t,
                       const PylithInt tindex,
                       const pylith::topology::Field& solution);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    pylith::string_vector _vertexDataFields; ///< Names of solution fields to output.

    static const char* _pyreComponent; ///< Name of Pyre component.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    OutputSoln(const OutputSoln&);   ///< Not implemented.
    const OutputSoln& operator=(const OutputSoln&);   ///< Not implemented

}; // OutputSoln

#endif // pylith_meshio_OutputSoln_hh

// End of file
