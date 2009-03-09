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

namespace pylith {
  namespace feassemble {

    template<typename quadrature_type>
    class Integrator
    { // Integrator

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Constructor
      Integrator(void);
      
      /// Destructor
      virtual
      ~Integrator(void);
      
      /** Set quadrature for integrating finite-element
       * quantities. Quadrature should already be initialized.
       *
       * @param q Quadrature for integrating.
       */
      void quadrature(const quadrature_type* q);
      
      /** Set manager of scales used to nondimensionalize problem.
       *
       * @param dim Nondimensionalizer.
       */
      void normalizer(const spatialdata::units::Nondimensional& dim);
      
      /** Set gravity field.
       *
       * @param g Gravity field.
       */
      void gravityField(spatialdata::spatialdb::GravityField* const gravityField);
      
      /** Set time step for advancing from time t to time t+dt.
       *
       * @param dt Time step
       */
      virtual
      void timeStep(const double dt);
      
      /** Get stable time step for advancing from time t to time t+dt.
       *
       * Default is MAXFLOAT (or 1.0e+30 if MAXFLOAT is not defined in math.h).
       *
       * @returns Time step
       */
      virtual
      double stableTimeStep(void) const;
      
      /** Check whether Jacobian needs to be recomputed.
       *
       * @returns True if Jacobian needs to be recomputed, false otherwise.
       */
      virtual
      bool needNewJacobian(void) const;
      
      /** Set flag for setting constraints for total field solution or
       *  incremental field solution.
       *
       * @param flag True if using incremental solution, false otherwise.
       */
      virtual
      void useSolnIncr(const bool flag);
      
      /** Initialize integrator.
       *
       * @param mesh Finite-element mesh.
       */
      virtual
      void initialize(const pylith::topology::Mesh& mesh);
      
      /** Integrate contributions to residual term (r) for operator.
       *
       * @param residual Field containing values for residual
       * @param t Current time
       * @param fields Solution fields
       */
      virtual 
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
      virtual
      void integrateJacobian(PetscMat* jacobian,
			     const double t,
			     pylith::topology::SolutionFields* const fields);

      /** Integrate contributions to residual term (r) for operator that
       * do not require assembly over cells, vertices, or processors.
       *
       * @param residual Field containing values for residual
       * @param t Current time
       * @param fields Solution fields
       */
      virtual 
      void integrateResidualAssembled(const pylith::topology::Field<pylith::topology::Mesh>& residual,
				      const double t,
				      pylith::topology::SolutionFields* const fields);

      /** Integrate contributions to Jacobian matrix (A) associated with
       * operator that do not require assembly over cells, vertices, or
       * processors
       *
       * @param jacobian Sparse matrix for Jacobian of system.
       * @param t Current time
       * @param fields Solution fields
       */
      virtual
      void integrateJacobianAssembled(PetscMat* jacobian,
				      const double t,
				      pylith::topology::SolutionFields* const fields);

      /** Update state variables as needed.
       *
       * @param t Current time
       * @param fields Solution fields
       * @param mesh Finite-element mesh
       */
      virtual
      void updateState(const double t,
		       pylith::topology::SolutionFields* const fields);

      /** Verify configuration is acceptable.
       *
       * @param mesh Finite-element mesh
       */
      virtual
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const = 0;


    }; // Integrator

  } // feassemble
} // pylith


// End of file 
