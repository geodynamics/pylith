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
