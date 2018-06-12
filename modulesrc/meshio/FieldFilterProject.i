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
 * @file modulesrc/meshio/FieldFilterProject.i
 *
 * @brief Python interface to C++ FieldFilterProject object.
 */

namespace pylith {
    namespace meshio {

        class pylith::meshio::FieldFilterProject : public FieldFilter {

	  // PUBLIC METHODS /////////////////////////////////////////////////
	public :
	  
	  /// Constructor
	  FieldFilterProject(void);
	  
	  /// Destructor
	  ~FieldFilterProject(void);
	  
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
	  
	}; // FieldFilterProject
      
    } // meshio
} // pylith


// End of file 
