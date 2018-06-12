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

        class pylith::feassemble::Observer : public pylith::utils::PyreComponent {
	  
	  // PUBLIC METHODS ///////////////////////////////////////////////////////
	public:
	  
	  /// Constructor
	  Observer(void);
	  
	  /// Destructor
	  virtual ~Observer(void);
	  
	  /** Check whether we want to write output at time t.
	   *
	   * @param[in] t Time of proposed write.
	   * @param[in] tindex Inxex of current time step.
	   * @returns True if output should be written at time t, false otherwise.
	   */
	  virtual
	  bool shouldWrite(const PylithReal t,
			   const PylithInt tindex) = 0;
	  
	}; // Observer

    } // feassemble
} // pylith


// End of file
