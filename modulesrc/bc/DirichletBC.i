// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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
		      const double upDir[3]);
      
      /** Set number of degrees of freedom that are constrained at
       * points in field.
       *
       * @param field Solution field
       */
      void setConstraintSizes(const pylith::topology::Field<pylith::topology::Mesh>& field);
      
      /** Set which degrees of freedom are constrained at points in field.
       *
       * @param field Solution field
       */
      void setConstraints(const pylith::topology::Field<pylith::topology::Mesh>& field);
      
      /** Set values in field.
       *
       * @param t Current time
       * @param field Solution field
       */
      void setField(const double t,
		    const pylith::topology::Field<pylith::topology::Mesh>& field);
      
      /** Set values in field.
       *
       * @param t0 Time t.
       * @param t1 Time t+dt.
       * @param field Solution field
       */
      void setFieldIncr(const double t0,
			const double t1,
			const pylith::topology::Field<pylith::topology::Mesh>& field);
      
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
