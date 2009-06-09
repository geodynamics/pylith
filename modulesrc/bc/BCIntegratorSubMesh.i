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

/** @file modulesrc/bc/BCIntegratorSubMesh.i
 *
 * @brief Python interface to C++ BCIntegratorSubMesh object.
 */

namespace pylith {
  namespace bc {

    class pylith::bc::BCIntegratorSubMesh : public BoundaryCondition,
		    public pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh> >
    { // class BoundaryCondition

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      BCIntegratorSubMesh(void);
      
      /// Destructor.
      virtual
      ~BCIntegratorSubMesh(void);
      
      /** Get parameter fields.
       *
       * @returns Parameter fields.
       */
      const pylith::topology::Fields<pylith::topology::Field<pylith::topology::Mesh> >*
      parameterFields(void) const;
      
      /** Get boundary mesh.
       *
       * @return Boundary mesh.
       */
      const pylith::topology::SubMesh& boundaryMesh(void) const;
      
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
