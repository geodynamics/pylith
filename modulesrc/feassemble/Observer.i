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
 * @file modulesrc/feassemble/Observer.i
 *
 * @brief Python interface to C++ Observer object.
 */

namespace pylith {
    namespace feassemble {

        class Observer {
	    
	    // PUBLIC METHODS ///////////////////////////////////////////////////////
	public:
	    
	    /// Constructor
	    Observer(void);
	    
	    /// Destructor
	    virtual ~Observer(void);
	    
	    /** Receive update (subject of observer).
	     *
	     * @param[in] t Current time.
	     * @param[in] tindex Current time step.
	     * @param[in] solution Solution at time t.
	     * @param[in] infoOnly Flag is true if this update is before solution is available (e.g., after initialization).
	     */
	    virtual
	    void update(const PylithReal t,
			const PylithInt tindex,
			const pylith::topology::Field& solution,
			const bool infoOnly=false) = 0;
	  
	}; // Observer

    } // feassemble
} // pylith


// End of file
