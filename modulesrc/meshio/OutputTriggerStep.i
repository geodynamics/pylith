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
 * @file modulesrc/meshio/OutputTriggerStep.i
 *
 * @brief Python interface to C++ OutputTriggerStep object.
 */

namespace pylith {
    namespace meshio {

        class pylith::meshio::OutputTriggerStep : public pylith::meshio::OutputTrigger {
	  
	  // PUBLIC METHODS ///////////////////////////////////////////////////////
	public:
	  
	  /// Constructor
	  OutputTriggerStep(void);
	  
	  /// Destructor
	  ~OutputTriggerStep(void);
	  
	  /** Check whether we want to write output at time t.
	   *
	   * @param[in] t Time of proposed write.
	   * @param[in] tindex Inxex of current time step.
	   * @returns True if output should be written at time t, false otherwise.
	   */
	  bool shouldWrite(const PylithReal t,
			   const PylithInt tindex);
	  
	  /** Set number of steps to skip between writes.
	   *
	   * @param[in] Number of steps to skip between writes.
	   */
	  void numStepsSkip(const int value);
	  
	  /** Get number of steps to skip between writes.
	   *
	   * @returns Number of steps to skip between writes.
	   */
	  int numStepsSkip(void) const;

	}; // OutputTriggerStep

    } // meshio
} // pylith


// End of file
