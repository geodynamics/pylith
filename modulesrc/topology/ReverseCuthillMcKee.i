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
