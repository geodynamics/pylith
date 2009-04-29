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

    class FaultCohesiveDyn : public FaultCohesive,
			     public pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh> >
    { // class FaultCohesiveDyn

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      FaultCohesiveDyn(void);
      
      /// Destructor.
      virtual
      ~FaultCohesiveDyn(void);
      
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

      /** Integrate contribution of cohesive cells to residual term.
       *
       * @param residual Field containing values for residual
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateResidual(const pylith::topology::Field<pylith::topology::Mesh>& residual,
			     const double t,
			     pylith::topology::SolutionFields* const fields);

      /** Integrate contributions to Jacobian matrix (A) associated with
       * operator.
       *
       * @param jacobian Sparse matrix for Jacobian of system.
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateJacobian(pylith::topology::Jacobian* jacobian,
			     const double t,
			     pylith::topology::SolutionFields* const fields);
  
      /** Verify configuration is acceptable.
       *
       * @param mesh Finite-element mesh
       */
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;
      
      /** Get vertex field associated with integrator.
       *
       * @param name Name of vertex field.
       * @param fields Solution fields.
       *
       * @returns Vertex field.
       */
      const pylith::topology::Field<pylith::topology::SubMesh>&
      vertexField(const char* name,
		  const pylith::topology::SolutionFields& fields);
      
      /** Get cell field associated with integrator.
       *
       * @param name Name of cell field.
       * @param fields Solution fields.
       *
       * @returns Cell field.
       */
      const pylith::topology::Field<pylith::topology::SubMesh>&
      cellField(const char* name,
		const pylith::topology::SolutionFields& fields);
      
    }; // class FaultCohesiveDyn

  } // faults
} // pylith


// End of file 
