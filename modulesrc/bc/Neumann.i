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

/** @file modulesrc/bc/Neumann.i
 *
 * @brief Python interface to C++ Neumann object.
 */

namespace pylith {
  namespace bc {

    class Neumann : public BCIntegratorSubMesh,
		    public TimeDependent
    { // class Neumann

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      Neumann(void);

      /// Destructor.
      ~Neumann(void);

      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
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
       * @param residual Field containing values for residual.
       * @param t Current time.
       * @param fields Solution fields.
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
      
      /** Get cell field with BC information.
       *
       * @param fieldType Type of field.
       * @param name Name of field.
       * @param mesh Finite-element mesh.
       * @param fields Solution fields.
       *
       * @returns Traction vector at integration points.
       */
      const pylith::topology::Field& cellField(const char* name,
					       pylith::topology::SolutionFields* const fields =0);

      // PROTECTED METHODS //////////////////////////////////////////////
    protected :
      
      /** Get label of boundary condition surface.
       *
       * @returns Label of surface (from mesh generator).
       */
      const char* _getLabel(void) const;
      
      /** Get manager of scales used to nondimensionalize problem.
       *
       * @returns Nondimensionalizer.
       */
      const spatialdata::units::Nondimensional& _getNormalizer(void) const;
      
    }; // class Neumann

  } // bc
} // pylith


// End of file 
