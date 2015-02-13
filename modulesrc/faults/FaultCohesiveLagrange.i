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

/** @file modulesrc/faults/FaultCohesiveLagrange.i
 *
 * @brief Python interface to C++ FaultCohesiveLagrange object.
 */

namespace pylith {
  namespace faults {

    class FaultCohesiveLagrange : public FaultCohesive
    { // class FaultCohesiveLagrange

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      FaultCohesiveLagrange(void);
      
      /// Destructor.
      virtual
      ~FaultCohesiveLagrange(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Initialize fault. Determine orientation and setup boundary
       * condition parameters.
       *
       * @param mesh Finite-element mesh.
       * @param upDir Direction perpendicular to along-strike direction that is 
       *   not collinear with fault normal (usually "up" direction but could 
       *   be up-dip direction; applies to fault surfaces in 2-D and 3-D).
       */
      virtual
      void initialize(const pylith::topology::Mesh& mesh,
		      const PylithScalar upDir[3]);
      
      /** Setup DOF on solution field.
       *
       * @param field Solution field.
       */
      void setupSolnDof(pylith::topology::Field* field);

      /** Integrate contributions to residual term (r) for operator that
       * do not require assembly across processors.
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
       * operator that do not require assembly processors.
       *
       * @param jacobian Sparse matrix
       * @param t Current time
       * @param fields Solution fields
       * @param mesh Finite-element mesh
       */
      void integrateJacobian(pylith::topology::Jacobian* jacobian,
			     const PylithScalar t,
			     pylith::topology::SolutionFields* const fields);
      
      /** Integrate contributions to Jacobian matrix (A) associated with
       * operator that do not require assembly across processors.
       *
       * @param jacobian Diagonal Jacobian matrix as a field.
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateJacobian(pylith::topology::Field* jacobian,
			     const PylithScalar t,
			     pylith::topology::SolutionFields* const fields);

      /** Adjust solution from solver with lumped Jacobian to match Lagrange
       *  multiplier constraints.
       *
       * @param fields Solution fields
       * @param t Current time.
       * @param jacobian Jacobian of the system.
       */
      virtual
      void adjustSolnLumped(pylith::topology::SolutionFields* fields,
			    const PylithScalar t,
			    const pylith::topology::Field& jacobian);

      /** Verify configuration is acceptable.
       *
       * @param mesh Finite-element mesh
       */
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;
      
      /** Verify constraints are acceptable.
       *
       * @param field Solution field.
       */
      virtual
      void checkConstraints(const pylith::topology::Field& solution) const;
      
    }; // class FaultCohesiveLagrange

  } // faults
} // pylith


// End of file 
