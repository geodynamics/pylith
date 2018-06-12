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
 * @file modulesrc/feassemble/ObservedSubject.i
 *
 * @brief Python interface to C++ ObservedSubject object.
 */

namespace pylith {
    namespace feassemble {

        class ObservedSubject {
	    
	    // PUBLIC METHODS ///////////////////////////////////////////////////////
	public:
	    /// Constructor.
	    ObservedSubject(void);
	    
	    /// Destructor
	    virtual ~ObservedSubject(void);
	    
	    /// Deallocate PETSc and local data structures.
	    virtual
	    void deallocate(void);
	    
	    /** Register observer to receive notifications.
	     *
	     * Observers are used for output.
	     *
	     * @param[in] observer Observer to receive notifications.
	     */
	    void registerObserver(pylith::feassemble::Observer* observer);
	    
	    /** Remove observer from receiving notifications.
	     *
	     * @param[in] observer Observer to remove.
	     */
	    void removeObserver(pylith::feassemble::Observer* observer);
	    
	    /** Send observers an update.
	     *
	     * @param[in] t Current time.
	     * @param[in] tindex Current time step.
	     * @param[in] solution Solution at time t.
	     * @param[in] infoOnly Flag is true if this update is before solution is available (e.g., after initialization).
	     */
	    void notifyObservers(const PylithReal t,
				 const PylithInt tindex,
				 const pylith::topology::Field& solution,
				 const bool infoOnly=false);
	    
	}; // ObservedSubject

    } // feassemble
} // pylith


// End of file
