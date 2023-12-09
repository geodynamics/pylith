// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================

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
