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
 * @file modulesrc/topology/ReverseCuthillMcKee.hh
 *
 * @brief Python interface to C++ PyLith ReverseCuthillMcKee object.
 */

namespace pylith {
  namespace topology {

    // ReverseCuthillMcKee ----------------------------------------------
    class ReverseCuthillMcKee
    { // ReverseCuthillMcKee

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /** Reorder vertices and cells of mesh using PETSc routines
       * implementing reverse Cuthill-McKee algorithm.
       *
       * @param mesh PyLith finite-element mesh.
       */
      static
      void reorder(topology::Mesh* mesh);

    }; // ReverseCuthillMcKee

  } // topology
} // pylith


// End of file
