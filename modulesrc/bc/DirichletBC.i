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
  namespace feassemble {
    class Constraint;
  } // feassemble

  namespace bc {

    class DirichletBC : public BoundaryCondition, 
			public pylith::feassemble::Constraint
    { // class DirichletBC

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      DirichletBC(void);
      
      /// Destructor.
      ~DirichletBC(void);
      
      /** Set database for rate of change of values.
       *
       * @param db Spatial database
       */
      void dbRate(spatialdata::spatialdb::SpatialDB* const db);
      
      /** Set indices of fixed degrees of freedom. 
       *
       * Note: all points associated with boundary condition has same
       * degrees of freedom fixed.
       *
       * Example: [0, 1] to fix x and y degrees of freedom in Cartesian system.
       *
       * @param flags Array of indices of fixed degrees of freedom.
       * @param size Size of array.
       */
      %apply(int* INPLACE_ARRAY1, int DIM1) {
	(const int* flags, 
	 const int size)
	  };
      void fixedDOF(const int* flags,
		    const int size);
      %clear(const int* flags, const int size);
      
      /** Set time at which rate of change begins.
       *
       * @param t Reference time.
       */
      void referenceTime(const double t);
      
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
      void setField(const double t,
		    const pylith::topology::Field& field);
      
    }; // class DirichletBC
    
  } // bc
} // pylith


// End of file 
