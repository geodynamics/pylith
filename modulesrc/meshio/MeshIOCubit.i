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
