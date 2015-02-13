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

      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Initialize boundary condition.
       *
       * @param mesh Finite-element mesh.
       * @param upDir Vertical direction (somtimes used in 3-D problems).
       */
      void initialize(const pylith::topology::Mesh& mesh,
		      const PylithScalar upDir[3]);

      /** Get boundary mesh.
       *
       * @return Boundary mesh.
       */
      const pylith::topology::Mesh& boundaryMesh(void) const;
      
      /** Get vertex field with BC information.
       *
       * @param name Name of field.
       * @param fields Solution fields.
       *
       * @returns Field over vertices.
       */
      const pylith::topology::Field&
      vertexField(const char* name,
		  const pylith::topology::SolutionFields& fields);
      
    }; // class DirichletBoundary
    
  } // bc
} // pylith


// End of file 
