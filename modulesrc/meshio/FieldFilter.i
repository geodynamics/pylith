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
 * @file modulesrc/meshio/FieldFilter.i
 *
 * @brief Python interface to C++ FieldFilter object.
 */

namespace pylith {
    namespace meshio {

        class pylith::meshio::FieldFilter : public pylith::utils::PyreComponent {

	  // PUBLIC METHODS ///////////////////////////////////////////////////////
	public:

	  /// Constructor
	  FieldFilter(void);
      
	  /// Destructor
	  virtual ~FieldFilter(void);
      
	  /// Deallocate PETSc and local data structures.
	  virtual
	  void deallocate(void);
      
	  /** Create copy of filter.
	   *
	   * @returns Copy of filter.
	   */
	  virtual
	  FieldFilter* clone(void) const = 0;
      
	  /** Filter field.
	   *
	   * @param fieldIn Field to filter.
	   * @returns Field after applying filter.
	   */
	  virtual
	  pylith::topology::Field* filter(pylith::topology::Field* fieldIn) = 0;
      
	}; // FieldFilter
    
    
    } // meshio
} // pylith


// End of file 
