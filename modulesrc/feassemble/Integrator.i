// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/feassemble/Integrator.i
 *
 * @brief Python interface to C++ abstract Integrator object.
 */

namespace pylith {
  namespace feassemble {

    class Integrator
    { // Integrator

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Constructor
      Integrator(void);
      
      /// Destructor
      virtual
      ~Integrator(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Get the quadrature for integrating finite-element
       * quantities.
       *
       * @returns Quadrature for integrating.
       */
      const pylith::feassemble::Quadrature& quadrature();
      
      /** Set quadrature for integrating finite-element
       * quantities. Quadrature should already be initialized.
       *
       * @param q Quadrature for integrating.
       */
      void quadrature(const pylith::feassemble::Quadrature* q);
      
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
      void timeStep(const PylithScalar dt);
      
      /** Get stable time step for advancing from time t to time t+dt.
       *
       * Default is MAXFLOAT (or 1.0e+30 if MAXFLOAT is not defined in math.h).
       *
       * @param mesh Finite-element mesh.
       * @returns Time step
       */
      virtual
      PylithScalar stableTimeStep(const pylith::topology::Mesh& mesh);
      
      /** Check whether Jacobian needs to be recomputed.
       *
       * @returns True if Jacobian needs to be recomputed, false otherwise.
       */
      virtual
      bool needNewJacobian(void) const;
      
      /** Check whether integrator generates a symmetric Jacobian.
       *
       * @returns True if integrator generates symmetric Jacobian.
       */
      virtual
      bool isJacobianSymmetric(void) const;

      /** Initialize integrator.
       *
       * @param mesh Finite-element mesh.
       */
      virtual
      void initialize(const pylith::topology::Mesh& mesh);
      
      /** Setup DOF on solution field.
       *
       * @param field Solution field.
       */
      virtual
      void setupSolnDof(pylith::topology::Field* field);

      /** Integrate contributions to residual term (r) for operator.
       *
       * @param residual Field containing values for residual
       * @param t Current time
       * @param fields Solution fields
       */
      virtual 
      void integrateResidual(const pylith::topology::Field& residual,
			     const PylithScalar t,
			     pylith::topology::SolutionFields* const fields);

      /** Integrate contributions to Jacobian matrix (A) associated with
       * operator.
       *
       * @param jacobian Sparse matrix for Jacobian of system.
       * @param t Current time
       * @param fields Solution fields
       */
      virtual
      void integrateJacobian(pylith::topology::Jacobian* jacobian,
			     const PylithScalar t,
			     pylith::topology::SolutionFields* const fields);

      /** Integrate contributions to Jacobian matrix (A) associated with
       * operator.
       *
       * @param jacobian Diagonal matrix (as field) for Jacobian of system.
       * @param t Current time
       * @param fields Solution fields
       */
      virtual
      void integrateJacobian(pylith::topology::Field* jacobian,
			     const PylithScalar t,
			     pylith::topology::SolutionFields* const fields);
      
      /** Update state variables as needed.
       *
       * @param t Current time
       * @param fields Solution fields
       * @param mesh Finite-element mesh
       */
      virtual
      void updateStateVars(const PylithScalar t,
			   pylith::topology::SolutionFields* const fields);

      /** Verify configuration is acceptable.
       *
       * @param mesh Finite-element mesh
       */
      virtual
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const = 0;

      /** Verify constraints are acceptable.
       *
       * @param field Solution field.
       */
      virtual
      void checkConstraints(const pylith::topology::Field& solution) const;


    }; // Integrator

  } // feassemble
} // pylith


// End of file 
