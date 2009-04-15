// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/**
 * @file modulesrc/problems/Formulation.hh
 *
 * @brief Python interface to C++ Formulation.
 */

namespace pylith {
  namespace problems {

    class Formulation
    { // Formulation
      
      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Constructor
      Formulation(void);

      /// Destructor
      ~Formulation(void);

      /** Set handles to integrators over the mesh.
       *
       * @param integrators Integrators over the mesh.
       * @param numIntegrators Number of integrators.
       */
      void meshIntegrators(pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::Mesh> >** integrators,
			   const int numIntegrators);
      
      /** Set handles to integrators over lower-dimension meshes.
       *
       * @param integrators Integrators over lower-dimension meshes.
       * @param numIntegrators Number of integrators.
       */
      void submeshIntegrators(pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh> >** integrators,
			      const int numIntegrators);
      
      /** Update handles and parameters for reforming the Jacobian and
       *  residual.
       *
       * @param jacobian Handle to sparse matrix for Jacobian of system.
       * @param fields Handle to solution fields.
       * @param t Current time (nondimensional).
       * @param dt Time step (nondimension).
       */
      void updateSettings(pylith::topology::Jacobian* jacobian,
			  pylith::topology::SolutionFields* fields,
			  const double t,
			  const double dt);
      
      /// Reform system residual.
      void reformResidual(PetscVec solutionVec, PetscVec residualVec);
      
      /// Reform system Jacobian.
      void reformJacobian(void);

    }; // Formulation

  } // problems
} // pylith


// End of file 
