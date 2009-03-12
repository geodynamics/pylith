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

/** @file modulesrc/bc/AbsorbingDampers.i
 *
 * @brief Python interface to C++ AbsorbingDampers object.
 */

%template(SubMeshIntegrator) pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh> >;

namespace pylith {
  namespace bc {

    class AbsorbingDampers : 
      public BoundaryCondition, 
      public pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh> >
    { // class AbsorbingDampers

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      AbsorbingDampers(void);

      /// Destructor.
      ~AbsorbingDampers(void);

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
       * @param mat Sparse matrix
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateJacobian(PetscMat* mat,
			     const double t,
			     pylith::topology::SolutionFields* const fields);
      
      /** Verify configuration is acceptable.
       *
       * @param mesh Finite-element mesh
       */
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;

    }; // class AbsorbingDampers

  } // bc
} //pylith


// End of file 
