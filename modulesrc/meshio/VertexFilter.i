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
 * @file modulesrc/meshio/VertexFilter.i
 *
 * @brief Python interface to C++ VertexFilter object.
 */

namespace pylith {
  namespace meshio {

    class pylith::meshio::VertexFilter
    { // VertexFilter

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      VertexFilter(void);
      
      /// Destructor
      virtual
      ~VertexFilter(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Create copy of filter.
       *
       * @returns Copy of filter.
       */
      virtual
      VertexFilter* clone(void) const = 0;
      
      /** Filter field over vertices of a mesh.
       *
       * @param fieldIn Field to filter.
       */
      virtual
      const pylith::topology::Field& filter(const pylith::topology::Field& fieldIn) = 0;

    }; // VertexFilter

  } // meshio
} // pylith


// End of file 
