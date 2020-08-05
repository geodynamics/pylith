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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/problems/TimeDependent.hh
 *
 * @brief Object for time dependent problem.
 */

#if !defined(pylith_problems_timedependent_hh)
#define pylith_problems_timedependent_hh

#include "Problem.hh" // ISA Problem
#include "pylith/testing/testingfwd.hh" // USES MMSTest

class pylith::problems::TimeDependent : public pylith::problems::Problem {
    friend class TestTimeDependent; // unit testing
    friend class pylith::testing::MMSTest; // Testing with Method of Manufactured Solutions

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TimeDependent(void);

    /// Destructor
    ~TimeDependent(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set start time for problem.
     *
     * @param[in] value Start time (seconds).
     */
    void setStartTime(const double value);

    /** Get start time for problem.
     *
     * @returns Start time (seconds).
     */
    double getStartTime(void) const;

    /** Set end time for problem.
     *
     * @param[in] value End time (seconds).
     */
    void setEndTime(const double value);

    /** Get end time for problem.
     *
     * @returns End time (seconds).
     */
    double getEndTime(void) const;

    /** Set maximum number of time steps.
     *
     * @param[in] value Maximum number of time steps.
     */
    void setMaxTimeSteps(const size_t value);

    /** Get maximum number of time steps.
     *
     * @returns Maximum number of time steps.
     */
    size_t getMaxTimeSteps(void) const;

    /** Set initial time step for problem.
     *
     * @param[in] value Initial time step (seconds).
     */
    void setInitialTimeStep(const double value);

    /** Get initial time step for problem.
     *
     * @returns Initial time step (seconds).
     */
    double getInitialTimeStep(void) const;

    /** Set initial conditions.
     *
     * @param[in] ic Array of initial conditions.
     * @param[in] numIC Number of initial conditions.
     */
    void setInitialCondition(pylith::problems::InitialCondition* ic[],
                             const int numIC);

    /** Should notify observers of solution with initial conditions.
     *
     * This will result in output being written at the starting time.
     *
     * @param[in] value True if observers should be notified of solution with initial conditions.
     */
    void setShouldNotifyIC(const bool value);

    /** Set progress monitor.
     *
     * @param[in] monitor Progress monitor for time-dependent simulation.
     */
    void setProgressMonitor(pylith::problems::ProgressMonitorTime* monitor);

    /** Get Petsc DM for problem.
     *
     * @returns PETSc DM for problem.
     */
    PetscDM getPetscDM(void);

    /** Get nonlinear solver for problem.
     *
     * @returns PETSc SNES for problem.
     */
    PetscSNES getPetscSNES(void);

    /// Verify configuration.
    void verifyConfiguration(void) const;

    /// Initialize.
    void initialize(void);

    /** Solve time dependent problem.
     */
    void solve(void);

    /** Perform operations before advancing solution one time step.
     *
     * Set constraints, etc.
     */
    void prestep(void);

    /** Perform Perform operations after advancing solution one time step.
     *
     * Update state variables, output.
     */
    void poststep(void);

    /** Callback static method for computing residual for RHS, G(t,s).
     *
     * @param[in] ts PETSc time stepper.
     * @param[in] t Current time.
     * @param[in] solutionVec PETSc Vec for solution.
     * @param[out] residualvec PETSc Vec for residual.
     * @param[in] context User context (TimeDependent).
     */
    static
    PetscErrorCode computeRHSResidual(PetscTS ts,
                                      PetscReal t,
                                      PetscVec solutionVec,
                                      PetscVec residualVec,
                                      void* context);

    /* Callback static method for computing Jacobian for RHS, Jacobian of G(t,s).
     *
     * @param[in] ts PETSc time stepper.
     * @param[in] t Current time.
     * @param[in] solution PetscVec for solution.
     * @param[out] jacobianMat Jacobian matrix.
     * @param[out] precondMat Preconditioner matrix.
     * @param[in] context User context (TimeDependent).
     */
    static
    PetscErrorCode computeRHSJacobian(PetscTS ts,
                                      PetscReal t,
                                      PetscVec solutionVec,
                                      PetscMat jacobianMat,
                                      PetscMat precondMat,
                                      void* context);

    /** Callback static method for computing residual for LHS, F(t,s,\dot{s}).
     *
     * @param[in] ts PETSc time stepper.
     * @param[in] t Current time.
     * @param[in] solutionVec PetscVec for solution.
     * @param[in] solutionDotVec PetscVec for time derivative of solution.
     * @param[out] residualvec PetscVec for residual.
     * @param[in] context User context (TimeDependent).
     */
    static
    PetscErrorCode computeLHSResidual(PetscTS ts,
                                      PetscReal t,
                                      PetscVec solutionVec,
                                      PetscVec solutionDotVec,
                                      PetscVec residualVec,
                                      void* context);

    /* Callback static method for computing Jacobian for LHS, Jacobian of F(t,s,\dot{s}).
     *
     * @param[in] ts PETSc time stepper.
     * @param[in] t Current time.
     * @param[in] solution PetscVec for solution.
     * @param[out] jacobianMat Jacobian matrix.
     * @param[out] precondMat Preconditioner matrix.
     * @param[in] context User context (TimeDependent).
     */
    static
    PetscErrorCode computeLHSJacobian(PetscTS ts,
                                      PetscReal t,
                                      PetscVec solutionVec,
                                      PetscVec solutionDotVec,
                                      PetscReal s_tshift,
                                      PetscMat jacobianMat,
                                      PetscMat precondMat,
                                      void* context);

    /** Callback static method for operations before advancing solution one time step.
     */
    static
    PetscErrorCode prestep(PetscTS ts);

    /** Callback static method for operations after advancing solution one time step.
     */
    static
    PetscErrorCode poststep(PetscTS ts);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    /// Notify observers with solution corresponding to initial conditions.
    void _notifyObserversInitialSoln(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    double _startTime; ///< Starting time of problem (seconds).
    double _endTime; ///< Ending time of problem (seconds).
    double _dtInitial; ///< Initial time step (seconds).
    size_t _maxTimeSteps; ///< Maximum number of time steps for problem.
    PetscTS _ts; ///< PETSc time stepper.
    std::vector<pylith::problems::InitialCondition*> _ic; ///< Array of initial conditions.
    pylith::problems::ProgressMonitorTime* _monitor; ///< Monitor for simulation progress.
    bool _shouldNotifyIC;

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    TimeDependent(const TimeDependent&); ///< Not implemented
    const TimeDependent& operator=(const TimeDependent&); ///< Not implemented

}; // TimeDependent

#endif // pylith_problems_timedependent_hh

// End of file
