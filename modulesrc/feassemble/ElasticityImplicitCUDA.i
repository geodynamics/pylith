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

/** @file modulesrc/feassemble/ElasticityImplicitCUDA.i
 *
 * @brief Python interface to C++ ElasticityImplicitCUDA object.
 */

namespace pylith {
  namespace feassemble {

    class ElasticityImplicitCUDA : public IntegratorElasticity
    { // ElasticityImplicitCUDA

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Constructor
      ElasticityImplicitCUDA(void);

      /// Destructor
      ~ElasticityImplicitCUDA(void);

      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Set time step for advancing from time t to time t+dt.
       *
       * @param dt Time step
       */
      void timeStep(const PylithScalar dt);
      
      /** Get stable time step for advancing from time t to time t+dt.
       *
       * Default is current time step.
       *
       * @param mesh Finite-element mesh.
       * @returns Time step
       */
      PylithScalar stableTimeStep(const pylith::topology::Mesh& mesh);
      
      /** Integrate residual part of RHS for 3-D finite elements.
       * Includes gravity and element internal force contribution.
       *
       * We assume that the effects of boundary conditions are already
       * included in the residual (tractions, concentrated nodal forces,
       * and contributions to internal force vector due to
       * displacement/velocity BC).  This routine computes the additional
       * external loads due to body forces plus the
       * element internal forces for the current stress state.
       *
       * @param residual Field containing values for residual
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateResidual(const pylith::topology::Field<pylith::topology::Mesh>& residual,
			     const PylithScalar t,
			     pylith::topology::SolutionFields* const fields);
      
      /** Integrate contributions to Jacobian matrix (A) associated with
       * operator.
       *
       * @param jacobian Sparse matrix for Jacobian of system.
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateJacobian(pylith::topology::Jacobian* jacobian,
			     const PylithScalar t,
			     pylith::topology::SolutionFields* const fields);
  
    }; // ElasticityImplicitCUDA

  } // feassemble
} // pylith


// End of file 
