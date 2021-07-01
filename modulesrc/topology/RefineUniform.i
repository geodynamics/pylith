// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/topology/RefineUniform.hh
 *
 * @brief Python interface to C++ PyLith RefineUniform object.
 */

namespace pylith {
  namespace topology {

    // RefineUniform ----------------------------------------------------
    class pylith::topology::RefineUniform
    { // RefineUniform

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Constructor
      RefineUniform(void);
      
      /// Destructor
      ~RefineUniform(void);
      
      /** Refine mesh.
       *
       * @param newMesh Refined mesh (result).
       * @param mesh Mesh to refine.
       * @param levels Number of levels to refine.
       */
      void refine(Mesh* const newMesh,
		  const Mesh& mesh,
		  const int levels =1);

    }; // RefineUniform

  } // topology
} // pylith


// End of file
