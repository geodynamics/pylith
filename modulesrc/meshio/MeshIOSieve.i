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
 * @file modulesrc/meshio/MeshIOSieve.i
 *
 * @brief Python interface to C++ MeshIOSieve object.
 */

namespace pylith {
  namespace meshio {

    class MeshIOSieve : public MeshIO
    { // MeshIOSieve

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      MeshIOSieve(void);

      /// Destructor
      ~MeshIOSieve(void);

      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Set filename for Sieve mesh file.
       *
       * @param filename Name of file
       */
      void filename(const char* name);
      
      /** Get filename of Sieve mesh file.
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

    }; // MeshIOSieve

  } // meshio
} // pylith


// End of file 
