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

/** @file modulesrc/feassemble/Integrator.i
 *
 * @brief Python interface to C++ abstract Integrator object.
 */

%template(MeshIntegrator) pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::Mesh> >;

namespace pylith {
  namespace feassemble {

    class IntegratorElasticity :
      public pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::Mesh> >
    { // IntegratorElasticity

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Constructor
      IntegratorElasticity(void);

      /// Destructor
      ~IntegratorElasticity(void);

      /** Set material.
       *
       * @param m Elastic material.
       */
      void material(pylith::materials::ElasticMaterial* m);
      
      /** Determine whether we need to recompute the Jacobian.
       *
       * @returns True if Jacobian needs to be recomputed, false otherwise.
       */
      bool needNewJacobian(void);
      
      /** Set flag for setting constraints for total field solution or
       *  incremental field solution.
       *
       * @param flag True if using incremental solution, false otherwise.
       */
      void useSolnIncr(const bool flag);
      
      /** Initialize integrator.
       *
       * @param mesh Finite-element mesh.
       */
      void initialize(const pylith::topology::Mesh& mesh);
      
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
      
    }; // IntegratorElasticity

  } // feassemble
} // pylith


// End of file 
