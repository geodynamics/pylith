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
 * @file libsrc/problems/GreensFns.hh
 *
 * @brief Object for Green's functions problem.
 */

#if !defined(pylith_problems_greensfns_hh)
#define pylith_problems_greensfns_hh

#define TEMPORARY

#include "Problem.hh" // ISA Problem
#include "pylith/testing/testingfwd.hh" // USES MMSTest
#include "pylith/faults/faultsfwd.hh" // HOLDSA FaultCohesiveImpulses

#include "pylith/testing/FaultCohesiveStub.hh" // TEMPORARY

class pylith::problems::GreensFns : public pylith::problems::Problem {
    friend class TestGreensFns; // unit testing
    friend class pylith::testing::MMSTest; // Testing with Method of Manufactured Solutions

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    GreensFns(void);

    /// Destructor
    ~GreensFns(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set id for fault with impulses.
     *
     * @param[in] value Id for fault with impulses.
     */
    void setFaultId(const int value);

    /** Get id for fault with impulses.
     *
     * @returns Id for fault with impulses.
     */
    int setFaultId(void) const;

    /** Set progress monitor.
     *
     * @param[in] monitor Progress monitor for Green's functions simulation.
     */
    void setProgressMonitor(pylith::problems::ProgressMonitorStep* monitor);

    /// Verify configuration.
    void verifyConfiguration(void) const;

    /// Initialize.
    void initialize(void);

    /** Solve Green's function problem.
     */
    void solve(void);

    /** Perform operations after advancing solution of one impulse
     *
     * Update state variables, output.
     */
    void poststep(void);

    /** Set solution values according to constraints (Dirichlet BC).
     *
     * @param[in] solutionVec PETSc Vec with current global view of solution.
     */
    void setSolutionLocal(PetscVec solutionVec);

    /** Compute residual and assemble into global vector.
     *
     * @param[out] residualVec PETSc Vec for residual.
     * @param[in] solutionVec PETSc Vec with current trial solution.
     */
    void computeResidual(PetscVec residualVec,
                         PetscVec solutionVec);

    /** Compute Jacobian for F(s,\dot{s})
     *
     * @param[out] jacobianMat PETSc Mat for Jacobian.
     * @param[out] precondMat PETSc Mat for preconditioner for Jacobian.
     * @param[in] solutionVec PETSc Vec with current trial solution.
     */
    void computeJacobian(PetscMat jacobianMat,
                         PetscMat precondMat,
                         PetscVec solutionVec);

    /** Callback static method for computing residual.
     *
     * @param[in] snes PETSc solver
     * @param[in] solutionVec PetscVec for solution.
     * @param[out] residualvec PetscVec for residual.
     * @param[in] context User context (GreensFns).
     */
    static
    PetscErrorCode computeResidual(PetscSNES snes,
                                      PetscVec solutionVec,
                                      PetscVec residualVec,
                                      void* context);

    /** Callback static method for computing Jacobian.
     *
     * @param[out] jacobianMat PETSc Mat for Jacobian.
     * @param[out] precondMat PETSc Mat for preconditioner for Jacobian.
     * @param[in] solutionVec PETSc Vec with current trial solution.
     */
    static
    void computeJacobian(PetscSNES snes,
                         PetscVec solutionVec,
                         PetscMat jacobianMat,
                         PetscMat precondMat,
                         void* context);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    PylithInt _faultImpulsesId;
#if defined(TEMPORARY)
    pylith::fa::FaultCohesiveStub _faultImpulses;
#else
    pylith::faults::FaultCohesiveImpulses* _faultImpulses; ///< Fault interface with Green's functions impulses.
#endif

    PetscSNES _snes; ///< PETSc SNES solver.
    pylith::problems::ProgressMonitorStep* _monitor; ///< Monitor for simulation progress.

    pylith::topology::Field* _residual; ///< Handle to residual field.
    pylith::topology::Field* _solutionDot; ///< Handle to time derivative of solution.

}; // GreensFns

#endif // pylith_problems_greensfns_hh

// End of file
