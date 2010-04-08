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

/** @file modulesrc/faults/FaultCohesiveLagrange.i
 *
 * @brief Python interface to C++ FaultCohesiveLagrange object.
 */

namespace pylith {
  namespace faults {

    class FaultCohesiveLagrange : public FaultCohesive
    { // class FaultCohesiveLagrange

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      FaultCohesiveLagrange(void);
      
      /// Destructor.
      virtual
      ~FaultCohesiveLagrange(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Initialize fault. Determine orientation and setup boundary
       * condition parameters.
       *
       * @param mesh Finite-element mesh.
       * @param upDir Direction perpendicular to along-strike direction that is 
       *   not collinear with fault normal (usually "up" direction but could 
       *   be up-dip direction; only applies to fault surfaces in a 3-D domain).
       * @param normalDir General preferred direction for fault normal
       *   (used to pick which of two possible normal directions for
       *   interface; only applies to fault surfaces in a 3-D domain).
       */
      virtual
      void initialize(const pylith::topology::Mesh& mesh,
		      const double upDir[3],
		      const double normalDir[3]);
      
      /** Split solution field for separate preconditioning.
       *
       * @param field Solution field.
       */
      void splitField(pylith::topology::Field<pylith::topology::Mesh>* field);

      /** Integrate contributions to residual term (r) for operator that
       * do not require assembly across cells, vertices, or processors.
       *
       * @param residual Field containing values for residual
       * @param t Current time
       * @param fields Solution fields
       */
      virtual
      void integrateResidualAssembled(const pylith::topology::Field<pylith::topology::Mesh>& residual,
				      const double t,
				      pylith::topology::SolutionFields* const fields);

      /** Integrate contributions to Jacobian matrix (A) associated with
       * operator that do not require assembly across cells, vertices, or
       * processors.
       *
       * @param jacobian Sparse matrix
       * @param t Current time
       * @param fields Solution fields
       * @param mesh Finite-element mesh
       */
      void integrateJacobianAssembled(pylith::topology::Jacobian* jacobian,
				      const double t,
				      pylith::topology::SolutionFields* const fields);
      
      /** Integrate contributions to Jacobian matrix (A) associated with
       * operator that do not require assembly across cells, vertices, or
       * processors.
       *
       * @param jacobian Diagonal Jacobian matrix as a field.
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateJacobianAssembled(pylith::topology::Field<pylith::topology::Mesh>* jacobian,
				      const double t,
				      pylith::topology::SolutionFields* const fields);

      /** Adjust solution from solver with lumped Jacobian to match Lagrange
       *  multiplier constraints.
       *
       * @param fields Solution fields
       * @param jacobian Jacobian of the system.
       */
      virtual
      void adjustSolnLumped(pylith::topology::SolutionFields* fields,
			    const pylith::topology::Field<pylith::topology::Mesh>& jacobian);

      /** Verify configuration is acceptable.
       *
       * @param mesh Finite-element mesh
       */
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;
      
      /** Verify constraints are acceptable.
       *
       * @param field Solution field.
       */
      virtual
      void checkConstraints(const pylith::topology::Field<pylith::topology::Mesh>& solution) const;
      
    }; // class FaultCohesiveLagrange

  } // faults
} // pylith


// End of file 
