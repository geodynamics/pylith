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

/** @file modulesrc/feassemble/ElasticityExplicit.i
 *
 * @brief Python interface to C++ ElasticityExplicit object.
 */

namespace pylith {
  namespace feassemble {

    class ElasticityExplicit : public IntegratorElasticity
    { // ElasticityExplicit

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :
      
      /// Constructor
      ElasticityExplicit(void);
      
      /// Destructor
      ~ElasticityExplicit(void);
      
      /** Set time step for advancing from time t to time t+dt.
       *
       * @param dt Time step
       */
      void timeStep(const double dt);
      
      /** Set flag for setting constraints for total field solution or
       *  incremental field solution.
       *
       * @param flag True if using incremental solution, false otherwise.
       */
      void useSolnIncr(const bool flag);
      
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
       * @param jacobian Sparse matrix for Jacobian of system.
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateJacobian(pylith::topology::Jacobian* jacobian,
			     const double t,
			     pylith::topology::SolutionFields* const fields);

    }; // ElasticityExplicit

  } // feassemble
} // pylith


// End of file 
