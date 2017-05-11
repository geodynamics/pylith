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
 * f(x,t) = f_0(x) + \dot{f}_1(x)*(t-t_1(x)) + f_2(x)*a(t-t_2(x)).
 *
 * Auxiliary fields:
 *     if _useInitial
 *         initial amplitude (scalar or vector) f_0(x)
 *    if _useRate
 *        rate amplitude (scalar or vector) \dot{f}_1(x)
 *        rate start (scalar) t_1(x)
 *    if _useTimeHistory
 *        time history amplitude (scalar or vector) f_2(x)
 *        time history start (scalar) t_2(x)
 *        time history value (scalar) a(t-t_2(x))
 */
class pylith::bc::DirichletTimeDependent : public DirichletNew { // class DirichletTimeDependent
    friend class TestDirichletTimeDependent;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    DirichletTimeDependent(void);

    /// Destructor.
    ~DirichletTimeDependent(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set time history database.
     *
     * @param[in] db Time history database.
     */
    void dbTimeHistory(spatialdata::spatialdb::TimeHistory* th);

    /** Get time history database.
     *
     * @preturns Time history database.
     */
    const spatialdata::spatialdb::TimeHistory* dbTimeHistory(void);

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

    /** Update auxiliary fields at beginning of time step.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     */
    void prestep(const double t,
                 const double dt);

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
     *
     * @param[in] solution Solution field.
     */
    void _auxFieldsSetup(const topology::Field& solution);

    /** Set kernels for RHS residual G(t,s).
     *
     * Potentially, there are g0 and g1 kernels for each equation. If no
     * kernel is needed, then set the kernel function to NULL.
     *
     * @param solution Solution field.
     */
    void _setFEKernelsConstraint(const topology::Field& solution);


    // PROTECTED MEMBERS //////////////////////////////////////////////////
private:

    spatialdata::spatialdb::TimeHistory* _dbTimeHistory;

    bool _useInitial; ///< Use initial value term.
    bool _useRate; ///< Use rate term.
    bool _useTimeHistory; ///< Use time history term.

    static const char* _pyreComponent; ///< Name of Pyre component.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    DirichletTimeDependent(const DirichletTimeDependent&); ///< Not implemented.
    const DirichletTimeDependent& operator=(const DirichletTimeDependent&); ///< Not implemented.

}; // class DirichletTimeDependent

#endif // pylith_bc_dirichlettimedependent_hh


// End of file
