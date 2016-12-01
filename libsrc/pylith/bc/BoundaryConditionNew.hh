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

#if !defined(pylith_bc_boundaryconditionnew_hh)
#define pylith_bc_boundaryconditionnew_hh

// Include directives ---------------------------------------------------
#include "bcfwd.hh" // forward declarations

#include <string> // HASA std::string

// BoundaryConditionNew ----------------------------------------------------
/// Abstract base class for boundary conditions.
class pylith::bc::BoundaryConditionNew
{ // class BoundaryConditionNew
friend class TestBoundaryConditionNew;   // unit testing

// PUBLIC METHODS /////////////////////////////////////////////////////
public:

/// Default constructor.
BoundaryConditionNew(void);

/// Destructor.
virtual
~BoundaryConditionNew(void);

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

// PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

std::string _label;   ///< Label to identify boundary condition points in mesh.

// NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

BoundaryConditionNew(const BoundaryConditionNew&); ///< Not implemented.
const BoundaryConditionNew& operator=(const BoundaryConditionNew&); ///< Not implemented.

}; // class BoundaryConditionNew

#endif // pylith_bc_boundaryconditionnew_hh


// End of file
