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
 * @file modulesrc/meshio/MeshIOAscii.i
 *
 * @brief Python interface to C++ MeshIOAscii object.
 */

namespace pylith {
  namespace meshio {

    class MeshIOAscii : public MeshIO
    { // MeshIOAscii

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      MeshIOAscii(void);

      /// Destructor
      ~MeshIOAscii(void);

      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Set filename for ASCII file.
       *
       * @param filename Name of file
       */
      void filename(const char* name);
      
      /** Get filename of ASCII file.
       *
       * @returns Name of file
       */
      const char* filename(void) const;

      // PROTECTED METHODS //////////////////////////////////////////////
    protected :

      /// Write mesh
      void _write(void) const;
      
      /// Read mesh
      void _read(void);

    }; // MeshIOAscii

  } // meshio
} // pylith


// End of file 
