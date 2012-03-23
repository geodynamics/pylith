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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
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
      
      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Set time step for advancing from time t to time t+dt.
       *
       * @param dt Time step
       */
      void timeStep(const PylithScalar dt);
      
      /** Set normalized viscosity for numerical damping.
       *
       * @param viscosity Nondimensional viscosity.
       */
      void normViscosity(const PylithScalar viscosity);

      /** Integrate contributions to residual term (r) for operator.
       *
       * @param residual Field containing values for residual
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateResidual(const pylith::topology::Field<pylith::topology::Mesh>& residual,
			     const PylithScalar t,
			     pylith::topology::SolutionFields* const fields);
      
      /** Integrate contributions to residual term (r) for operator.
       *
       * @param residual Field containing values for residual
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateResidualLumped(const pylith::topology::Field<pylith::topology::Mesh>& residual,
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

      /** Integrate contributions to Jacobian matrix (A) associated
       * with operator that require assembly across cells, vertices,
       * or processors.
       *
       * @param jacobian Diagonal Jacobian matrix as a field.
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateJacobian(pylith::topology::Field<pylith::topology::Mesh>* jacobian,
			     const PylithScalar t,
			     pylith::topology::SolutionFields* const fields);

    }; // ElasticityExplicit

  } // feassemble
} // pylith


// End of file 
