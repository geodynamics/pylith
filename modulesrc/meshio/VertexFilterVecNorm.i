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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/meshio/VertexFilterVecNorm.i
 *
 * @brief Python interface to C++ VertexFilterVecNorm object.
 */

namespace pylith {
  namespace meshio {

    template<typename field_type>
    class pylith::meshio::VertexFilterVecNorm : public VertexFilter<field_type>
    { // VertexFilterVecNorm

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      VertexFilterVecNorm(void);

      /// Destructor
      ~VertexFilterVecNorm(void);
      
      /** Create copy of filter.
       *
       * @returns Copy of filter.
       */
      VertexFilter<field_type>* clone(void) const;
      
      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Filter vertex field.
       *
       * @param fieldIn Field to filter.
       */
      const field_type& filter(const field_type& fieldIn);
      
    }; // VertexFilterVecNorm

  } // meshio
} // pylith


// End of file 
