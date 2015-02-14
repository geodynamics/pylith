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

/** @file modulesrc/feassemble/Constraint.i
 *
 * @brief Python interface to C++ abstract base Constraint.
 */

namespace pylith {
  namespace feassemble {

    class Constraint
    { // class Constraint

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      Constraint(void);

      /// Destructor.
      virtual
      ~Constraint(void);

      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Get number of constraints per location.
       *
       * @returns Number of constraints per location.
       */
      virtual
      int numDimConstrained(void) const = 0;

      /** Set manager of scales used to nondimensionalize problem.
       *
       * @param dim Nondimensionalizer.
       */
      void normalizer(const spatialdata::units::Nondimensional& dim);

      /** Set number of degrees of freedom that are constrained at
       * points in field.
       *
       * @param field Solution field
       */
      virtual
      void setConstraintSizes(const pylith::topology::Field& field) = 0;

      /** Set which degrees of freedom are constrained at points in field.
       *
       * @param field Solution field
       */
      virtual
      void setConstraints(const pylith::topology::Field& field) = 0;

      /** Set values in field.
       *
       * @param t Current time
       * @param field Solution field
       */
      virtual
      void setField(const PylithScalar t,
		    const pylith::topology::Field& field) = 0;
      
      /** Set increment in values from t0 to t1 in field.
       *
       * @param t0 Time t.
       * @param t1 Time t+dt.
       * @param field Solution field
       */
      virtual
      void setFieldIncr(const PylithScalar t0,
			const PylithScalar t1,
			const pylith::topology::Field& field) = 0;

    }; // class Constraint

  } // feassemble
} // pylith


// End of file 
