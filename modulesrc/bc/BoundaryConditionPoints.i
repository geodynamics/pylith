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

/** @file modulesrc/bc/BoundaryConditionPoints.i
 *
 * @brief Python interface to C++ BoundaryConditionPoints object.
 */

namespace pylith {
  namespace bc {

    class pylith::bc::BoundaryConditionPoints : public BoundaryCondition
    { // class BoundaryCondition

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      BoundaryConditionPoints(void);
      
      /// Destructor.
      virtual
      ~BoundaryConditionPoints(void);

      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Get parameter fields.
       *
       * @returns Parameter fields.
       */
      const pylith::topology::Fields* parameterFields(void) const;

    }; // class BoundaryConditionPoints

  } // bc
} // pylith


// End of file 
