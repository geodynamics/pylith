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

/** @file libsrc/bc/Neumann.hh
 *
 * @brief C++ implementation of Neumann (e.g., traction) boundary conditions.
 */

#if !defined(pylith_bc_neumann_hh)
#define pylith_bc_neumann_hh

// Include directives ---------------------------------------------------
#include "pylith/feassemble/IntegratorBoundary.hh" // ISA IntegratorBoundary

#include "pylith/bc/bcfwd.hh" // forward delcarations

#include "pylith/topology/topologyfwd.hh" // USES Field

// Neumann ----------------------------------------------------
/// @brief Neumann (e.g., traction) boundary conditions.
class pylith::bc::Neumann : public pylith::feassemble::IntegratorBoundary {
    friend class TestNeumann;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    Neumann(void);

    /// Destructor.
    ~Neumann(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Name of scale associated with Neumann boundary
     * condition (e.g., 'pressure' for elasticity).
     *
     * A Neumann boundary condition constrains the gradient in
     * a solution subfield. In some cases the constraint is
     * actually on a scaled version of the gradient as is the
     * case of a Neumann boundary condition for elasticity
     * that constrains boundary tractions.
     *
     * @param value Name of scale for nondimensionalizing Neumann boundary condition.
     */
    void scaleName(const char* value);


    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::FieldBase::Description _description; ///< Description of field associated with BC.
    std::string _scaleName; ///< Name of scale associated with Neumann boundary condition.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    Neumann(const Neumann&); ///< Not implemented.
    const Neumann& operator=(const Neumann&); ///< Not implemented.

}; // class Neumann

#endif // pylith_bc_neumann_hh


// End of file
