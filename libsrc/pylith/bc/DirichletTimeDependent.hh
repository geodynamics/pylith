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
// See LICENSE.md.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/bc/DirichletTimeDependent.hh
 *
 * @brief C++ implementation of Dirichlet (prescribed values at
 * degrees of freedom) boundary condition with time-dependent expression
 * for value.
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

#if !defined(pylith_bc_dirichlettimedependent_hh)
#define pylith_bc_dirichlettimedependent_hh

#include "pylith/bc/BoundaryCondition.hh" // ISA BoundaryCondition

#include "pylith/topology/topologyfwd.hh" // USES Field

class pylith::bc::DirichletTimeDependent : public pylith::bc::BoundaryCondition {
    friend class TestDirichletTimeDependent; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    DirichletTimeDependent(void);

    /// Destructor.
    ~DirichletTimeDependent(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set indices of constrained degrees of freedom at each location.
     *
     * Example: [0, 1] to apply forces to x and y degrees of freedom in
     * a Cartesian coordinate system.
     *
     * @param[in] dof Array of indices for constrained degrees of freedom.
     * @param[in] size Size of array
     */
    void setConstrainedDOF(const int* flags,
                           const int size);

    /** Get indices of constrained degrees of freedom.
     *
     * @returns Array of indices for constrained degrees of freedom.
     */
    const pylith::int_array& getConstrainedDOF(void) const;

    /** Set time history database.
     *
     * @param[in] db Time history database.
     */
    void setTimeHistoryDB(spatialdata::spatialdb::TimeHistory* th);

    /** Get time history database.
     *
     * @preturns Time history database.
     */
    const spatialdata::spatialdb::TimeHistory* getTimeHistoryDB(void);

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

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Integrator if applicable, otherwise NULL.
     */
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution);

    /** Create constraint and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Constraint if applicable, otherwise NULL.
     */
    pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution);

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& domainMesh);

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Derived field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& domainMesh);

    /** Update time-dependent auxiliary field.
     *
     * @param[out] auxiliaryField Auxiliary field.
     * @param[in] t Current time.
     */
    void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                              const double t);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    /** Update kernel constants.
     *
     * @param[in] dt Current time step.
     */
    void _updateKernelConstants(const PylithReal dt);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    int_array _constrainedDOF; ///< List of constrained degrees of freedom at each location.
    spatialdata::spatialdb::TimeHistory* _dbTimeHistory; ///< Time history database.
    pylith::bc::TimeDependentAuxiliaryFactory* _auxiliaryFactory; ///< Factory for auxiliary subfields.

    bool _useInitial; ///< Use initial value term.
    bool _useRate; ///< Use rate term.
    bool _useTimeHistory; ///< Use time history term.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    DirichletTimeDependent(const DirichletTimeDependent&); ///< Not implemented.
    const DirichletTimeDependent& operator=(const DirichletTimeDependent&); ///< Not implemented.

}; // class DirichletTimeDependent

#endif // pylith_bc_dirichlettimedependent_hh

// End of file
