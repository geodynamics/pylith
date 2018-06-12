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
 * @file libsrc/problems/ProblemObserver.hh
 *
 * @brief Observer of problem.
 */

#if !defined(pylith_problems_problemobserver_hh)
#define pylith_problems_problemobserver_hh

// Include directives ---------------------------------------------------
#include "problemsfwd.hh" // forward declarations

#include "pylith/feassemble/Observer.hh" // ISA Observer

// ProblemObserver --------------------------------------------------------
/// Observer of integrator for governing equations. Receives updates of solution.
class pylith::problems::ProblemObserver : public pylith::feassemble::Observer {
    friend class TestProblemObserver;   // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /** Constructor with problem to observe.
     *
     * Holds problem and does not manage its memory.
     *
     * @param[in] problem Problem to observe.
     */
    ProblemObserver(pylith::problems::Problem* const problem);


    /// Destructor
    virtual ~ProblemObserver(void);


    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);


    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    pylith::problems::Problem* const _problem;

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    ProblemObserver(const ProblemObserver&);   ///< Not implemented.
    const ProblemObserver& operator=(const ProblemObserver&);   ///< Not implemented

}; // ProblemObserver

#endif // pylith_problems_problemobserver_hh


// End of file
