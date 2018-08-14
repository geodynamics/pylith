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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/meshio/FieldFilterNone.i
 *
 * @brief Python interface to C++ FieldFilterNone object.
 */

namespace pylith {
    namespace meshio {

        class pylith::meshio::FieldFilterNone : public FieldFilter {

	  // PUBLIC METHODS /////////////////////////////////////////////////
	public :
	  
	  /// Constructor
	  FieldFilterNone(void);
	  
	  /// Destructor
	  ~FieldFilterNone(void);
	  
	  /** Create copy of filter.
	   *
	   * @returns Copy of filter.
	   */
	  FieldFilter* clone(void) const;
	  
	  /// Deallocate PETSc and local data structures.
	  void deallocate(void);
	  
	  /** Filter vertex field.
	   *
	   * @param fieldIn Field to filter.
	   */
	  const pylith::topology::Field* filter(pylith::topology::Field* fieldIn);
	  
	}; // FieldFilterNone
      
    } // meshio
} // pylith


// End of file 
