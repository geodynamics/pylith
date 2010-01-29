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

/** @file modulesrc/feassemble/ElasticityExplicitLgDeform.i
 *
 * @brief Python interface to C++ ElasticityExplicitLgDeform object.
 */

namespace pylith {
  namespace feassemble {

    class ElasticityExplicitLgDeform : public IntegratorElasticityLgDeform
    { // ElasticityExplicitLgDeform

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public:

      /// Constructor
      ElasticityExplicitLgDeform(void);

      /// Destructor
      ~ElasticityExplicitLgDeform(void);

      /// Deallocate PETSc and local data structures.
      void deallocate(void);
      
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

      /** Integrate contributions to Jacobian matrix (A) associated
       * with operator that require assembly across cells, vertices,
       * or processors.
       *
       * @param jacobian Diagonal Jacobian matrix as a field.
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateJacobian(pylith::topology::Field<pylith::topology::Mesh>* jacobian,
			     const double t,
			     pylith::topology::SolutionFields* const fields);

    }; // ElasticityExplicitLgDeform

  } // feassemble
} // pylith
