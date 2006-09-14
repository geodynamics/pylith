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

#if !defined(pylith_meshio_meshio_hh)
#define pylith_meshio_meshio_hh

#include <Mesh.hh> // PETSc Mesh

namespace pylith {
  namespace meshIO {
    class MeshIO;
  } // meshio
} // pylith

class pylith::meshIO::MeshIO
{ // MeshIO
  
// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  MeshIO(void);

  /// Destructor
  virtual ~MeshIO(void);

  /** Read mesh from file.
   *
   * @param pMesh Pointer to PETSc mesh object
   */
  virtual void read(ALE::Obj<ALE::Mesh>& pMesh) = 0;

  /** Write mesh to file.
   *
   * @param mesh PETSc mesh object
   */
  virtual void write(const ALE::Obj<ALE::Mesh>& mesh) const = 0;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  /** Get flag indicating whether indices start at 0 (True) or 1 (False).
   *
   * @returns True if indices start at 0, false if 1.
   */
  bool useIndexZero(void) const;

  /** Set flag indicating whether indices start at 0 (True) or 1 (False).
   *
   * @param flag True if indices start at 0, false if 1.
   */
  void useIndexZero(const bool flag);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  bool _useIndexZero; ///< Flag indicating if indicates start at 0 (T) or 1 (F)

}; // MeshIO

#include "MeshIO.icc" // inline methods

#endif // pylith_meshio_meshio_hh

// End of file 
