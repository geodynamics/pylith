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

/** @file modulesrc/bc/AbsorbingDampers.i
 *
 * @brief Python interface to C++ AbsorbingDampers object.
 */

namespace pylith {
  namespace bc {

    class AbsorbingDampers : public BCIntegratorSubMesh
    { // class AbsorbingDampers

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      AbsorbingDampers(void);

      /// Destructor.
      ~AbsorbingDampers(void);

      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Set database for boundary condition parameters.
       *
       * @param db Spatial database
       */
      void db(spatialdata::spatialdb::SpatialDB* const db);

      /** Initialize boundary condition.
       *
       * @param mesh Finite-element mesh.
       * @param upDir Direction perpendicular to horizontal surface tangent 
       *   direction that is not collinear with surface normal.
       */
      void initialize(const pylith::topology::Mesh& mesh,
		      const PylithScalar upDir[3]);
      
      /** Integrate contributions to residual term (r) for operator.
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
       * @param mat Sparse matrix
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateJacobian(pylith::topology::Jacobian* mat,
			     const PylithScalar t,
			     pylith::topology::SolutionFields* const fields);
      
      /** Integrate contributions to Jacobian matrix (A) associated with
       * operator.
       *
       * @param jacobian Jacobian of system.
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateJacobian(pylith::topology::Field* jacobian,
			     const PylithScalar t,
			     pylith::topology::SolutionFields* const fields);

      /** Verify configuration is acceptable.
       *
       * @param mesh Finite-element mesh
       */
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;

    }; // class AbsorbingDampers

  } // bc
} //pylith


// End of file 
