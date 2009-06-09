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

/** @file modulesrc/bc/Neumann.i
 *
 * @brief Python interface to C++ Neumann object.
 */

namespace pylith {
  namespace bc {

    class Neumann : public BCIntegratorSubMesh
    { // class Neumann

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      Neumann(void);

      /// Destructor.
      ~Neumann(void);

      /** Set database for boundary condition parameters.
       *
       * @param db Spatial database
       */
      void db(spatialdata::spatialdb::SpatialDB* const db);

      /** Initialize boundary condition.
       *
       * @param mesh Finite-element mesh.
       * @param upDir Direction perpendicular to horizontal surface tangent 
       *   direction that is not collinear with surface normal.
       */
      void initialize(const pylith::topology::Mesh& mesh,
		      const double upDir[3]);

      /** Integrate contributions to residual term (r) for operator.
       *
       * @param residual Field containing values for residual.
       * @param t Current time.
       * @param fields Solution fields.
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
      
      /** Get cell field with BC information.
       *
       * @param fieldType Type of field.
       * @param name Name of field.
       * @param mesh Finite-element mesh.
       * @param fields Solution fields.
       *
       * @returns Traction vector at integration points.
       */
      const pylith::topology::Field<pylith::topology::SubMesh>&
      cellField(const char* name,
		pylith::topology::SolutionFields* const fields =0);

    }; // class Neumann

  } // bc
} // pylith


// End of file 
