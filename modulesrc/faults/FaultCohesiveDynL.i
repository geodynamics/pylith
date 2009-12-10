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

/** @file modulesrc/faults/FaultCohesiveDynL.i
 *
 * @brief Python interface to C++ FaultCohesiveDynL object.
 */

namespace pylith {
  namespace faults {

    class FaultCohesiveDynL : public FaultCohesive
    { // class FaultCohesiveDynL

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      FaultCohesiveDynL(void);
      
      /// Destructor.
      virtual
      ~FaultCohesiveDynL(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
      
      /** Sets the spatial database for the inital tractions.
       *
       * @param db spatial database for initial tractions
       */
      void dbInitial(spatialdata::spatialdb::SpatialDB* db);
  
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
      void initialize(const pylith::topology::Mesh& mesh,
		      const double upDir[3],
		      const double normalDir[3]);
      
      /** Split solution field for separate preconditioning.
       *
       * @param field Solution field.
       */
      void splitField(pylith::topology::Field<pylith::topology::Mesh>* field);

      /** Integrate contributions to residual term (r) for operator that
       * require assembly across processors.
       *
       * @param residual Field containing values for residual
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateResidual(const pylith::topology::Field<pylith::topology::Mesh>& residual,
			     const double t,
			     pylith::topology::SolutionFields* const fields);

      /** Integrate contributions to residual term (r) for operator that
       * do not require assembly across cells, vertices, or processors.
       *
       * @param residual Field containing values for residual
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateResidualAssembled(const pylith::topology::Field<pylith::topology::Mesh>& residual,
				      const double t,
				      pylith::topology::SolutionFields* const fields);

      /** Integrate contributions to Jacobian matrix (A) associated with
       * operator that do not require assembly across cells, vertices, or
       * processors.
       *
       * @param mat Sparse matrix
       * @param t Current time
       * @param fields Solution fields
       * @param mesh Finite-element mesh
       */
      void integrateJacobianAssembled(pylith::topology::Jacobian* jacobian,
				      const double t,
				      pylith::topology::SolutionFields* const fields);
      
      /** Update state variables as needed.
       *
       * @param t Current time
       * @param fields Solution fields
       * @param mesh Finite-element mesh
       */
      void updateStateVars(const double t,
			   pylith::topology::SolutionFields* const fields);
      
      /** Constrain solution space based on friction.
       *
       * @param fields Solution fields.
       * @param t Current time.
       * @param jacobian Sparse matrix for system Jacobian.
       */
      void constrainSolnSpace(pylith::topology::SolutionFields* const fields,
			      const double t,
			      const pylith::topology::Jacobian& jacobian);

      /** Verify configuration is acceptable.
       *
       * @param mesh Finite-element mesh
       */
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;
      
      /** Get vertex field associated with integrator.
       *
       * @param name Name of cell field.
       * @param fields Solution fields.
       * @returns Vertex field.
       */
      const pylith::topology::Field<pylith::topology::SubMesh>&
      vertexField(const char* name,
		  const pylith::topology::SolutionFields* fields =0);
      
      /** Get cell field associated with integrator.
       *
       * @param name Name of cell field.
       * @param fields Solution fields.
       * @returns Cell field.
       */
      const pylith::topology::Field<pylith::topology::SubMesh>&
      cellField(const char* name,
		const pylith::topology::SolutionFields* fields =0);

      /** Cohesive cells use Lagrange multiplier constraints?
       *
       * @returns True if implementation using Lagrange multiplier
       * constraints, false otherwise.
       */
      bool useLagrangeConstraints(void) const;

      /** Get fields associated with fault.
       *
       * @returns Fields associated with fault.
       */
      const pylith::topology::Fields<pylith::topology::Field<pylith::topology::SubMesh> >*
      fields(void) const;

    }; // class FaultCohesiveDynL

  } // faults
} // pylith


// End of file 
