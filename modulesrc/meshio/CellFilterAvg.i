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
      
      /** Filter field over cells.
       *
       * @param fieldIn Field to filter.
       * @param label Label identifying cells.
       * @param labelId Value of label of cells to filter.
       *
       * @returns Averaged field.
       */
      const field_type& filter(const field_type& fieldIn,
			       const char* label =0,
			       const int labelId =0);

    }; // CellFilterAvg

  } // meshio
} // pylith


// End of file 
