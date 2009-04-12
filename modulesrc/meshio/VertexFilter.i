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
 * @file modulesrc/meshio/VertexFilter.i
 *
 * @brief Python interface to C++ VertexFilter object.
 */

namespace pylith {
  namespace meshio {

    template<typename field_type>
    class pylith::meshio::VertexFilter
    { // VertexFilter

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      VertexFilter(void);
      
      /// Destructor
      ~VertexFilter(void);
      
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
      const field_type& filter(const field_type& fieldIn) = 0;

    }; // VertexFilter

  } // meshio
} // pylith


// End of file 
