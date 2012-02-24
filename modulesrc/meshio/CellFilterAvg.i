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
 * @file modulesrc/meshio/CellFilterAvg.i
 *
 * @brief Python interface to C++ CellFilterAvg object.
 */

namespace pylith {
  namespace meshio {

    template<typename mesh_type, typename field_type>
    class pylith::meshio::CellFilterAvg :
      public CellFilter<mesh_type, field_type>
    { // CellFilterAvg

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      CellFilterAvg(void);

      /// Destructor
      ~CellFilterAvg(void);

      /** Create copy of filter.
       *
       * @returns Copy of filter.
       */
      CellFilter<mesh_type, field_type>* clone(void) const;
      
      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Get averaged field buffer.
       *
       * @returns Field associated with averaged values.
       */
      const field_type* fieldAvg(void) const;
  
      /** Filter field over cells.
       *
       * @param fieldIn Field to filter.
       * @param label Label identifying cells.
       * @param labelId Value of label of cells to filter.
       *
       * @returns Averaged field.
       */
      field_type& filter(const field_type& fieldIn,
			 const char* label =0,
			 const int labelId =0);

    }; // CellFilterAvg

  } // meshio
} // pylith


// End of file 
