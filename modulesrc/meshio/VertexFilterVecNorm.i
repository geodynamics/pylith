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
 * @file modulesrc/meshio/VertexFilterVecNorm.i
 *
 * @brief Python interface to C++ VertexFilterVecNorm object.
 */

namespace pylith {
  namespace meshio {

    template<typename mesh_type>
    class pylith::meshio::VertexFilterVecNorm : public VertexFilter<mesh_type>
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
      VertexFilter<mesh_type>* clone(void) const;
      
      /** Filter vertex field.
       *
       * @param fieldIn Field to filter.
       */
      const pylith::topology::Field<mesh_type>&
      filter(const pylith::topology::Field<mesh_type>& fieldIn);
      
    }; // VertexFilterVecNorm

  } // meshio
} // pylith


// End of file 
