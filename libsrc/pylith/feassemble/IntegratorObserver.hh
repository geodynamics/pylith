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
 * @file libsrc/feassemble/IntegratorObserver.hh
 *
 * @brief Observer of integrator for governing equations.
 */

#if !defined(pylith_feassemble_integratorobserver_hh)
#define pylith_feassemble_integratorobserver_hh

// Include directives ---------------------------------------------------
#include "feassemblefwd.hh" // forward declarations

#include "pylith/feassemble/Observer.hh" // ISA Observer

// IntegratorObserver --------------------------------------------------------
/// Observer of integrator for governing equations. Receives updates of solution.
class pylith::feassemble::IntegratorObserver : public pylith::feassemble::Observer {
    friend class TestIntegratorObserver;   // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /** Constructor with integrator to observe.
     *
     * Holds integrator and does not manage its memory.
     *
     * @param[in] integrator Integrator to observe.
     */
    IntegratorObserver(pylith::feassemble::IntegratorPointwise* const integrator);


    /// Destructor
    virtual ~IntegratorObserver(void);


    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);


    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    pylith::feassemble::IntegratorPointwise* const _integrator;

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    IntegratorObserver(const IntegratorObserver&);   ///< Not implemented.
    const IntegratorObserver& operator=(const IntegratorObserver&);   ///< Not implemented

}; // IntegratorObserver

#endif // pylith_feassemble_integratorobserver_hh


// End of file
