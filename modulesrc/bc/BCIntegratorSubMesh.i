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

/** @file modulesrc/bc/BCIntegratorSubMesh.i
 *
 * @brief Python interface to C++ BCIntegratorSubMesh object.
 */

namespace pylith {
  namespace bc {

    class pylith::bc::BCIntegratorSubMesh : public BoundaryCondition,
					    public pylith::feassemble::Integrator
    { // class BoundaryCondition

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      BCIntegratorSubMesh(void);
      
      /// Destructor.
      virtual
      ~BCIntegratorSubMesh(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Get parameter fields.
       *
       * @returns Parameter fields.
       */
      const pylith::topology::Fields* parameterFields(void) const;
      
      /** Get boundary mesh.
       *
       * @return Boundary mesh.
       */
      const pylith::topology::Mesh& boundaryMesh(void) const;
      
      /** Get mesh labels for submesh associated with applied forces.
       *
       * @param mesh Finite-element mesh.
       */
      void createSubMesh(const pylith::topology::Mesh& mesh);
      
      /** Verify configuration is acceptable.
       *
       * @param mesh Finite-element mesh
       */
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;
      
    }; // BCIntegratorSubMesh

  } // bc
} // pylith


// End of file 
