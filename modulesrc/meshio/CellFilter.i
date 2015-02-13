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
 * @file modulesrc/meshio/CellFilter.i
 *
 * @brief Python interface to C++ CellFilter object.
 */

namespace pylith {
  namespace meshio {

    class pylith::meshio::CellFilter
    { // CellFilter

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :
      
      /// Constructor
      CellFilter(void);
      
      /// Destructor
      virtual
      ~CellFilter(void);
      
      /** Create copy of filter.
       *
       * @returns Copy of filter.
       */
      virtual
      CellFilter* clone(void) const = 0;
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set quadrature associated with cells.
       *
       * @param q Quadrature for cells.
       */
      void quadrature(const pylith::feassemble::Quadrature* q);

      /** Filter field. Field type of filtered field is returned via an argument.
       *
       * @param fieldIn Field to filter.
       * @param label Value of label of cells to filter.
       * @param labelId Id associated with label of cells to filter.
       *
       * @returns Averaged field.
       */
      virtual
      pylith::topology::Field& filter(const pylith::topology::Field& fieldIn,
				      const char* label =0,
				      const int labelId =0) = 0;

    }; // CellFilter

  } // meshio
} // pylith


// End of file 
