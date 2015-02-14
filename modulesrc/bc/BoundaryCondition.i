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

/** @file modulesrc/bc/BoundaryCondition.i
 *
 * @brief Python interface to C++ BoundaryCondition object.
 */

namespace pylith {
  namespace bc {

    class BoundaryCondition
    { // class BoundaryCondition

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      BoundaryCondition(void);

      /// Destructor.
      virtual
      ~BoundaryCondition(void);

      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set label of boundary condition surface.
       *
       * @param value Label of surface (from mesh generator).
       */
      void label(const char* value);

      /** Get label of boundary condition surface.
       *
       * @returns Label of surface (from mesh generator).
       */
      const char* label(void) const;

      /** Verify configuration.
       *
       * @param mesh Finite-element mesh.
       */
      virtual
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;

      /** Initialize boundary condition.
       *
       * @param mesh Finite-element mesh.
       * @param upDir Vertical direction (somtimes used in 3-D problems).
       */
      virtual
      void initialize(const pylith::topology::Mesh& mesh,
		      const PylithScalar upDir[3]) = 0;

    }; // class BoundaryCondition

  } // bc
} // pylith


// End of file 
