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

/** @file modulesrc/bc/DirichletBC.i
 *
 * @brief Python interface to C++ DirichletBC object.
 */

namespace pylith {
  namespace bc {

    class DirichletBC : public TimeDependentPoints, 
			public pylith::feassemble::Constraint
    { // class DirichletBC

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      DirichletBC(void);
      
      /// Destructor.
      virtual
      ~DirichletBC(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Get number of constraints per location.
       *
       * @returns Number of constraints per location.
       */
      int numDimConstrained(void) const;

      /** Initialize boundary condition.
       *
       * @param mesh PETSc mesh
       * @param upDir Vertical direction (somtimes used in 3-D problems).
       */
      void initialize(const pylith::topology::Mesh& mesh,
		      const PylithScalar upDir[3]);
      
      /** Set number of degrees of freedom that are constrained at
       * points in field.
       *
       * @param field Solution field
       */
      void setConstraintSizes(const pylith::topology::Field& field);
      
      /** Set which degrees of freedom are constrained at points in field.
       *
       * @param field Solution field
       */
      void setConstraints(const pylith::topology::Field& field);
      
      /** Set values in field.
       *
       * @param t Current time
       * @param field Solution field
       */
      void setField(const PylithScalar t,
		    const pylith::topology::Field& field);
      
      /** Set values in field.
       *
       * @param t0 Time t.
       * @param t1 Time t+dt.
       * @param field Solution field
       */
      void setFieldIncr(const PylithScalar t0,
			const PylithScalar t1,
			const pylith::topology::Field& field);
      
      // PROTECTED METHODS //////////////////////////////////////////////
    protected :
      
      /** Get manager of scales used to nondimensionalize problem.
       *
       * @returns Nondimensionalizer.
       */
      const spatialdata::units::Nondimensional& _getNormalizer(void) const;

    }; // class DirichletBC
    
  } // bc
} // pylith


// End of file 
