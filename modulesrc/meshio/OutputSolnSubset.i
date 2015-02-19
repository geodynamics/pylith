// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/**
 * @file modulesrc/meshio/OutputSolnSubset.i
 *
 * @brief Python interface to C++ OutputSolnSubset object.
 */

namespace pylith {
  namespace meshio {

    class pylith::meshio::OutputSolnSubset : public OutputManager
    { // OutputSolnSubset

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      OutputSolnSubset(void);
      
      /// Destructor
      ~OutputSolnSubset(void);
      
      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
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
      const pylith::topology::Mesh& subdomainMesh(const pylith::topology::Mesh& mesh);
  
    }; // OutputSolnSubset

  } // meshio
} // pylith


// End of file 
