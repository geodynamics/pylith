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
 * @file modulesrc/meshio/OutputTriggerTime.i
 *
 * @brief Python interface to C++ OutputTriggerTime object.
 */

namespace pylith {
    namespace meshio {

        class pylith::meshio::OutputTriggerTime : public pylith::meshio::OutputTrigger {

	  // PUBLIC METHODS ///////////////////////////////////////////////////////
	public:

	  /// Constructor
	  OutputTriggerTime(void);

	  /// Destructor
	  ~OutputTriggerTime(void);

	  /** Check whether we want to write output at time t.
	   *
	   * @param[in] t Time of proposed write.
	   * @param[in] tindex Inxex of current time step.
	   * @returns True if output should be written at time t, false otherwise.
	   */
	  bool shouldWrite(const PylithReal t,
			   const PylithInt tindex);

	  /** Set elapsed time between writes.
	   *
	   * @param[in] Elapsed time between writes.
	   */
	  void setTimeSkip(const double value);

	  /** Get elapsed time between writes.
	   *
	   * @returns Elapsed time between writes.
	   */
	  double getTimeSkip(void) const;
	}; // OutputTriggerTime

    } // meshio
} // pylith


// End of file
