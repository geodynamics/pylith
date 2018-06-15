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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/OutputSolnBoundary.hh
 *
 * @brief C++ object for managing output of finite-element data over a
 * boundary.
 */

#if !defined(pylith_meshio_outputsolnboundary_hh)
#define pylith_meshio_outputsolnboundary_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/meshio/OutputSoln.hh" // ISA OutputSoln

#include "pylith/topology/topologyfwd.hh" // HOLDSA Mesh

#include <string> // HASA std::string

// OutputSolnBoundary -----------------------------------------------------
/** @brief C++ object for managing output of finite-element data over
 * a boundary.
 */
class pylith::meshio::OutputSolnBoundary : public pylith::meshio::OutputSoln {

    friend class TestOutputSolnBoundary; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /** Constructor
     *
     * @param[in] problem Problem to observe.
     */
    OutputSolnBoundary(pylith::problems::Problem* const problem);

    /// Destructor
    virtual ~OutputSolnBoundary(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set label identifier for subdomain.
     *
     * @param[in] value Label of subdomain.
     */
    void label(const char* value);

    /** Verify configuration.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    // PROTECTED METHODS ////////////////////////////////////////////////////
protected:

    /** Write solution at time step.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     */
    virtual
    void _writeDataStep(const PylithReal t,
                        const PylithInt tindex,
                        const pylith::topology::Field& solution);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    std::string _label; ///< Label of subdomain.
    pylith::topology::Mesh* _boundaryMesh; ///< Mesh of subdomain.

    static const char* _pyreComponent; ///< Name of Pyre component.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    OutputSolnBoundary(const OutputSolnBoundary&); ///< Not implemented.
    const OutputSolnBoundary& operator=(const OutputSolnBoundary&); ///< Not implemented

}; // OutputSolnBoundary

#endif // pylith_meshio_outputsolnboundary_hh

// End of file
