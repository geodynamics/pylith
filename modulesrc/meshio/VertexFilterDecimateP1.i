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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/meshio/VertexFilterDecimateP1.i
 *
 * @brief Python interface to C++ VertexFilterDecimateP1 object.
 */

namespace pylith {
  namespace meshio {

    class pylith::meshio::VertexFilterDecimateP1 : public VertexFilter
    { // VertexFilterDecimateP1

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      VertexFilterDecimateP1(void);

      /// Destructor
      ~VertexFilterDecimateP1(void);
      
      /** Create copy of filter.
       *
       * @returns Copy of filter.
       */
      VertexFilter* clone(void) const;
      
      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Filter vertex field.
       *
       * @param fieldIn Field to filter.
       */
      const pylith::topology::Field& filter(const pylith::topology::Field& fieldIn);
      
    }; // VertexFilterDecimateP1

  } // meshio
} // pylith


// End of file 
