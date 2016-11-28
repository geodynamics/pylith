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

/** @file libsrc/bc/DirichletTimeDependent.hh
 *
 * @brief C++ implementation of Dirichlet (prescribed values at
 * degrees of freedom) boundary condition with time-dependent expression
 * for value.
 */

#if !defined(pylith_bc_dirichlettimedependent_hh)
#define pylith_bc_dirichlettimedependent_hh

// Include directives ---------------------------------------------------
#include "DirichletNew.hh" // ISA DirichletNew

// DirichletTimeDependent ----------------------------------------------------
/** @brief Dirichlet (prescribed values at degrees of freedom) boundary
 * condition with a time-depdnent expression.
 *
 * f(x,t) = f_0(x) + \dot{f}_0(x)*(t-t_0(x)) + f_1(x)*a(t-t_1(x)).
 */
class pylith::bc::DirichletTimeDependent : public DirichletNew
{ // class DirichletTimeDependent
friend class TestDirichletTimeDependent;   // unit testing

// PUBLIC METHODS /////////////////////////////////////////////////////
public:

/// Default constructor.
DirichletTimeDependent(void);

/// Destructor.
~DirichletTimeDependent(void);

/// Deallocate PETSc and local data structures.
void deallocate(void);

/** Use initial value term in time history expression.
 *
 * @param[in] value True if using initial value term in expression.
 */
void useInitial(const bool value);

/** Get flag associated with using initial value term in time history expression.
 *
 * @returns True if using initial value term in expression, false otherwise.
 */
bool useInitial(void) const;

/** Use rate value term in time history expression.
 *
 * @param[in] value True if using rate value term in expression.
 */
void useRate(const bool value);

/** Get flag associated with using rate value term in time history expression.
 *
 * @returns True if using rate value term in expression, false otherwise.
 */
bool useRate(void) const;

/** Use time history term in time history expression.
 *
 * @param[in] value True if using time history term in expression.
 */
void useTimeHistory(const bool value);

/** Get flag associated with using time history term in time history expression.
 *
 * @returns True if using time history term in expression, false otherwise.
 */
bool useTimeHistory(void) const;

// PROTECTED METHODS //////////////////////////////////////////////////
protected:

/** Setup auxiliary subfields (discretization and query fns).
 *
 * Create subfields in auxiliary fields (includes name of the field,
 * vector field type, discretization, and scale for
 * nondimensionalization) and set query functions for filling them
 * from a spatial database.
 *
 * @attention The order of the calls to subfieldAdd() must match the
 * order of the auxiliary fields in the FE kernels.
 */
void _auxFieldsSetup(void);

/** Set kernels for RHS residual G(t,s).
 *
 * Potentially, there are g0 and g1 kernels for each equation. If no
 * kernel is needed, then set the kernel function to NULL.
 *
 * @param solution Solution field.
 */
void _setFEKernelsConstraint(const topology::Field& solution) const;


// PROTECTED MEMBERS //////////////////////////////////////////////////
private:

bool _useInitial;
bool _useRate;
bool _useTimeHistory;

// NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

DirichletTimeDependent(const DirichletTimeDependent&); ///< Not implemented.
const DirichletTimeDependent& operator=(const DirichletTimeDependent&); ///< Not implemented.

}; // class DirichletTimeDependent

#endif // pylith_bc_dirichlettimedependent_hh


// End of file
