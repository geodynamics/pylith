// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

/** @file libsrc/bc/BoundaryCondition.hh
 *
 * @brief C++ abstract base class for BoundaryCondition object.
 *
 * Interface definition for boundary conditions.
 */

#if !defined(pylith_bc_boundarycondition_hh)
#define pylith_bc_boundarycondition_hh

// Include directives ---------------------------------------------------
#include "bcfwd.hh" // forward declarations

#include <string> // HASA std::string

// BoundaryCondition ----------------------------------------------------
/// Abstract base class for boundary conditions.
class pylith::bc::BoundaryCondition { // class BoundaryCondition
    friend class TestBoundaryCondition;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    BoundaryCondition(void);

    /// Destructor.
    virtual ~BoundaryCondition(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set mesh label associated with boundary condition surface.
     *
     * @param[in] value Label of surface (from mesh generator).
     */
    void label(const char* value);

    /** Get mesh label associated with boundary condition surface.
     *
     * @returns Label of surface (from mesh generator).
     */
    const char* label(void) const;

    /** Set name of field in solution to constrain.
     *
     * @param[in] value Name of field in solution to constrain.
     */
    void field(const char* value);

    /** Get name of field in solution to constrain.
     *
     * @returns Name of field in solution to constrain.
     */
    const char* field(void) const;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    std::string _label;   ///< Label to identify boundary condition points in mesh.
    std::string _field; ///< Name of solution field for boundary condition.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    BoundaryCondition(const BoundaryCondition&); ///< Not implemented.
    const BoundaryCondition& operator=(const BoundaryCondition&); ///< Not implemented.

}; // class BoundaryCondition

#endif // pylith_bc_boundarycondition_hh


// End of file
