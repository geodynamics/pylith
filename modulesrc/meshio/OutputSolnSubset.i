// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/**
 * @file modulesrc/meshio/OutputSolnSubset.i
 *
 * @brief Python interface to C++ OutputSolnSubset object.
 */

%template(_SubMeshOutputManager) pylith::meshio::OutputManager<pylith::topology::SubMesh>;

namespace pylith {
  namespace meshio {

    class pylith::meshio::OutputSolnSubset : public OutputManager<pylith::topology::SubMesh>
    { // OutputSolnSubset

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      OutputSolnSubset(void);
      
      /// Destructor
      ~OutputSolnSubset(void);
      
      /** Set label identifier for subdomain.
       *
       * @param value Label of subdomain.
       */
      void label(const char* value);
      
      /** Verify configuration.
       *
       * @param mesh PETSc mesh
       */
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;

      /** Get mesh associated with subdomain.
       *
       * @returns Mesh associated with subdomain.
       */
      const pylith::topology::SubMesh& subdomainMesh(const pylith::topology::Mesh& mesh);
  
      /** Append finite-element vertex field to file.
       *
       * @param t Time associated with field.
       * @param field Vertex field.
       */
      void appendVertexField(const double t,
			     const pylith::topology::Field<pylith::topology::Mesh>& field);



    }; // OutputSolnSubset

  } // meshio
} // pylith


// End of file 
