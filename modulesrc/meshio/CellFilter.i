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
 * @file modulesrc/meshio/CellFilter.i
 *
 * @brief Python interface to C++ CellFilter object.
 */

namespace pylith {
  namespace meshio {

    template<typename mesh_type, typename field_type>
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
      cirtual
      void deallocate(void);
  
      /** Set quadrature associated with cells.
       *
       * @param q Quadrature for cells.
       */
      void quadrature(const pylith::feassemble::Quadrature<mesh_type>* q);

      /** Filter field. Field type of filtered field is returned via an argument.
       *
       * @param fieldIn Field to filter.
       * @param label Value of label of cells to filter.
       * @param labelId Id associated with label of cells to filter.
       *
       * @returns Averaged field.
       */
      virtual
      const field_type& filter(const field_type& fieldIn,
			       const char* label =0,
			       const int labelId =0) = 0;

    }; // CellFilter

  } // meshio
} // pylith


// End of file 
