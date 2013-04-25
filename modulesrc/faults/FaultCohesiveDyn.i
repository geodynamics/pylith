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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/faults/FaultCohesiveDyn.i
 *
 * @brief Python interface to C++ FaultCohesiveDyn object.
 */

namespace pylith {
  namespace faults {

    class FaultCohesiveDyn : public FaultCohesiveLagrange
    { // class FaultCohesiveDyn

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      FaultCohesiveDyn(void);
      
      /// Destructor.
      virtual
      ~FaultCohesiveDyn(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
      
      /** Sets the traction perturbation for prescribed tractions.
       *
       * @param tract Spatial and temporal variation of tractions.
       */
      void tractPerturbation(TractPerturbation* tract);
  
      /** Set the friction (constitutive) model.
       *
       * @param model Fault constutive model.
       */
      void frictionModel(pylith::friction::FrictionModel* const model);

      /** Nondimensional tolerance for detecting near zero values.
       *
       * @param value Nondimensional tolerance
       */
      void zeroTolerance(const PylithScalar value);

      /** Set flag used to determine when fault is traction free when it
       * opens or it still imposes any initial tractions.
       *
       * If true, acts as a frictional contact. If false, one can simulate
       * a dike opening.
       *
       * @param value Nondimensional tolerance
       */
      void openFreeSurf(const bool value);

      /** Initialize fault. Determine orientation and setup boundary
       * condition parameters.
       *
       * @param mesh Finite-element mesh.
       * @param upDir Direction perpendicular to along-strike direction that is 
       *   not collinear with fault normal (usually "up" direction but could 
       *   be up-dip direction; applies to fault surfaces in 2-D and 3-D).
       */
      void initialize(const pylith::topology::Mesh& mesh,
		      const PylithScalar upDir[3]);
      
      /** Integrate contributions to residual term (r) for operator that
       * do not require assembly across processors.
       *
       * Initial tractions (if specified) are already assembled and
       * contribute to the residual like Neumann boundary conditions.
       *
       * @param residual Field containing values for residual
       * @param t Current time
       * @param fields Solution fields
       */
      virtual
      void integrateResidual(const pylith::topology::Field<pylith::topology::Mesh>& residual,
				      const PylithScalar t,
				      pylith::topology::SolutionFields* const fields);

      /** Update state variables as needed.
       *
       * @param t Current time
       * @param fields Solution fields
       * @param mesh Finite-element mesh
       */
      void updateStateVars(const PylithScalar t,
			   pylith::topology::SolutionFields* const fields);
      
      /** Constrain solution space based on friction.
       *
       * @param fields Solution fields.
       * @param t Current time.
       * @param jacobian Sparse matrix for system Jacobian.
       */
      void constrainSolnSpace(pylith::topology::SolutionFields* const fields,
			      const PylithScalar t,
			      const pylith::topology::Jacobian& jacobian);

      /** Adjust solution from solver with lumped Jacobian to match Lagrange
       *  multiplier constraints.
       *
       * @param fields Solution fields.
       * @param t Current time.
       * @param jacobian Jacobian of the system.
       */
      void adjustSolnLumped(pylith::topology::SolutionFields* fields,
			    const PylithScalar t,
			    const pylith::topology::Field<pylith::topology::Mesh>& jacobian);

      /** Verify configuration is acceptable.
       *
       * @param mesh Finite-element mesh
       */
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;
      
      /** Get vertex field associated with integrator.
       *
       * @param name Name of cell field.
       * @param fields Solution fields.
       * @returns Vertex field.
       */
      const pylith::topology::Field<pylith::topology::SubMesh>&
      vertexField(const char* name,
		  const pylith::topology::SolutionFields* fields =0);
      
      /** Get cell field associated with integrator.
       *
       * @param name Name of cell field.
       * @param fields Solution fields.
       * @returns Cell field.
       */
      const pylith::topology::Field<pylith::topology::SubMesh>&
      cellField(const char* name,
		const pylith::topology::SolutionFields* fields =0);

    }; // class FaultCohesiveDyn

  } // faults
} // pylith


// End of file 
