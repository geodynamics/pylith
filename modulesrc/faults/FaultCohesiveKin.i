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

/** @file modulesrc/faults/FaultCohesiveKin.i
 *
 * @brief Python interface to C++ FaultCohesiveKin object.
 */

namespace pylith {
  namespace faults {

    class FaultCohesiveKin : public FaultCohesive,
			     public pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh> >
    { // class FaultCohesiveKin

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      FaultCohesiveKin(void);
      
      /// Destructor.
      virtual
      ~FaultCohesiveKin(void);
      
      /** Set kinematic earthquake sources.
       *
       * @param names Array of kinematic earthquake source names.
       * @param sources Array of kinematic earthquake sources.
       * @param numSources Number of earthquake sources
       */
      void eqsrcs(const char** names,
		  EqKinSrc** sources,
		  const int numSources);
      
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
       * @param matDB Database of bulk elastic properties for fault region
       *   (used to improve conditioning of Jacobian matrix)
       */
      void initialize(const pylith::topology::Mesh& mesh,
		      const double upDir[3],
		      const double normalDir[3],
		      spatialdata::spatialdb::SpatialDB* matDB);
      
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
		  const pylith::topology::SolutionFields& fields);
      
      /** Get cell field associated with integrator.
       *
       * @param name Name of cell field.
       * @param fields Solution fields.
       * @returns Cell field.
       */
      const pylith::topology::Field<pylith::topology::SubMesh>&
      cellField(const char* name,
		const pylith::topology::SolutionFields& fields);

    }; // class FaultCohesiveKin

  } // faults
} // pylith


// End of file 
