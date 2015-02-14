// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/bc/TimeDependentPoints.i
 *
 * @brief Python interface to C++ TimeDependentPoints object.
 */

namespace pylith {
  namespace bc {

    class pylith::bc::TimeDependentPoints : public BoundaryConditionPoints, 
					    public TimeDependent
    { // class TimeDependentPoints

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :
      
      /// Default constructor.
      TimeDependentPoints(void);
      
      /// Destructor.
      ~TimeDependentPoints(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set indices of degrees of freedom associated with BC.
       *
       * Note: Forces at all points are applied to the same degrees of freedom.
       *
       * Example: [0, 1] to apply forces to x and y degrees of freedom in
       * Cartesian system.
       *
       * @param flags Array of indices for degrees of freedom for forces.
       * @param size Size of array
       */
      %apply(int* INPLACE_ARRAY1, int DIM1) {
	(const int* flags, 
	 const int size)
	  };
      void bcDOF(const int* flags,
		 const int size);  
      %clear(const int* flags, const int size);
      
      // PROTECTED METHODS //////////////////////////////////////////////
    protected :
      
      /** Get label of boundary condition surface.
       *
       * @returns Label of surface (from mesh generator).
       */
      const char* _getLabel(void) const;

    }; // class TimeDependentPoints

  } // bc
} // pylith


// End of file 
