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

/** @file modulesrc/faults/FaultCohesiveTract.i
 *
 * @brief Python interface to C++ FaultCohesiveTract object.
 */

namespace pylith {
  namespace faults {

    class FaultCohesiveTract : public FaultCohesive
    { // class FaultCohesiveTract

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      FaultCohesiveTract(void);
      
      /// Destructor.
      virtual
      ~FaultCohesiveTract(void);
      
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
      void initialize(const pylith::topology::Mesh& mesh,
		      const PylithScalar upDir[3]);

      /** Integrate contribution of cohesive cells to residual term.
       *
       * @param residual Field containing values for residual
       * @param t Current time
       * @param fields Solution fields
       */
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
      void integrateJacobian(pylith::topology::Jacobian* jacobian,
			     const PylithScalar t,
			     pylith::topology::SolutionFields* const fields);
  
      /** Verify configuration is acceptable.
       *
       * @param mesh Finite-element mesh
       */
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;
      
      /** Get vertex field associated with integrator.
       *
       * @param name Name of vertex field.
       * @param fields Solution fields.
       *
       * @returns Vertex field.
       */
      const pylith::topology::Field& vertexField(const char* name,
						 const pylith::topology::SolutionFields* fields =0);
      
      /** Get cell field associated with integrator.
       *
       * @param name Name of cell field.
       * @param fields Solution fields.
       *
       * @returns Cell field.
       */
      const pylith::topology::Field& cellField(const char* name,
					       const pylith::topology::SolutionFields* fields =0);

    }; // class FaultCohesiveTract

  } // faults
} // pylith


// End of file 
