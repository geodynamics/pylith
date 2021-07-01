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
 * @file modulesrc/meshio/MeshIOCubit.i
 *
 * @brief Python interface to C++ MeshIOCubit object.
 */

namespace pylith {
  namespace meshio {

    class MeshIOCubit : public MeshIO
    { // MeshIOCubit

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      MeshIOCubit(void);

      /// Destructor
      ~MeshIOCubit(void);

      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Set filename for Cubit file.
       *
       * @param filename Name of file
       */
      void filename(const char* name);
      
      /** Get filename of Cubit file.
       *
       * @returns Name of file
       */
      const char* filename(void) const;
      
      /** Set flag on whether to use nodeset ids or names.
       *
       * @param flag True to use node set names.
       */
      void useNodesetNames(const bool flag);

      // PROTECTED METHODS ////////////////////////////////////////////////////
    protected :
      
      /// Write mesh
      void _write(void) const;
      
      /// Read mesh
      void _read(void);
      
    }; // MeshIOCubit

  } // meshio
} // pylith


// End of file 
