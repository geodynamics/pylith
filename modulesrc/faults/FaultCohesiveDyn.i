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

/** @file modulesrc/faults/FaultCohesiveDyn.i
 *
 * @brief Python interface to C++ FaultCohesiveDyn object.
 */

namespace pylith {
  namespace faults {

    class FaultCohesiveDyn : public FaultCohesiveLagrange
    { // class FaultCohesiveDyn

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      FaultCohesiveDyn(void);
      
      /// Destructor.
      virtual
      ~FaultCohesiveDyn(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
      
      /** Sets the spatial database for the inital tractions.
       *
       * @param db spatial database for initial tractions
       */
      void dbInitialTract(spatialdata::spatialdb::SpatialDB* db);
  
      /** Set the friction (constitutive) model.
       *
       * @param model Fault constutive model.
       */
      void frictionModel(pylith::friction::FrictionModel* const model);

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
      
      /** Integrate contributions to residual term (r) for operator that
       * do not require assembly across processors.
       *
       * Initial tractions (if specified) are already assembled and
       * contribute to the residual like Neumann boundary conditions.
       *
       * @param residual Field containing values for residual
       * @param t Current time
       * @param fields Solution fields
       */
      virtual
      void integrateResidual(const pylith::topology::Field<pylith::topology::Mesh>& residual,
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

      /** Adjust solution from solver with lumped Jacobian to match Lagrange
       *  multiplier constraints.
       *
       * @param fields Solution fields.
       * @param jacobian Jacobian of the system.
       */
      void adjustSolnLumped(pylith::topology::SolutionFields* fields,
			    const pylith::topology::Field<pylith::topology::Mesh>& jacobian);

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

    }; // class FaultCohesiveDyn

  } // faults
} // pylith


// End of file 
