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

/** @file modulesrc/bc/DirichletBoundary.i
 *
 * @brief Python interface to C++ DirichletBoundary object.
 */

namespace pylith {
  namespace bc {

    class DirichletBoundary : public DirichletBC
    { // class DirichletBoundary

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      DirichletBoundary(void);

      /// Destructor.
      ~DirichletBoundary(void);

      /** Initialize boundary condition.
       *
       * @param mesh Finite-element mesh.
       * @param upDir Vertical direction (somtimes used in 3-D problems).
       */
      void initialize(const pylith::topology::Mesh& mesh,
		      const double upDir[3]);

      /** Get boundary mesh.
       *
       * @return Boundary mesh.
       */
      const ALE::Obj<pylith::SieveSubMesh>& boundaryMesh(void) const;
      
      /** Get vertex field with BC information.
       *
       * @param fieldType Type of field.
       * @param name Name of field.
       * @param mesh Finite-element mesh.
       * @param fields Solution fields.
       *
       * @returns Field over vertices.
       */
      const pylith::topology::Field&
      vertexField(const char* name,
		  const pylith::topology::Mesh& mesh,
		  const pylith::topology::SolutionFields& fields);
      
    }; // class DirichletBoundary
    
  } // bc
} // pylith


// End of file 
