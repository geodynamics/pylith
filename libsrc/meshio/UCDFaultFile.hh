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
 * @file libsrc/meshio/UCDFaultFile.hh
 *
 * @brief C++ object for reading a fault mesh from a UCD file.
 *
 * Temporary fix for getting a fault mesh information into
 * PyLith. Using a fault mesh permits cells to have more than one face
 * on the fault surface. Supporting information also provides the
 * orientation of the fault surface, eliminating the need to determine
 * it from the faces/vertices alone.
 */

#if !defined(pylith_meshio_ucdfaultfile_hh)
#define pylith_meshio_ucdfaultfile_hh

// Include directives ---------------------------------------------------
#define NEWPYLITHMESH 1 
#include "pylith/utils/sievetypes.hh" // USES SieveMesh, ALE::Mesh

// Forward declarations -------------------------------------------------
namespace pylith {
  namespace meshio {
    class UCDFaultFile;

    class TestFaultFile; // unit testing
  } // meshio
} // pylith

// MeshIOLagrit ---------------------------------------------------------
class pylith::meshio::UCDFaultFile
{ // UCDFaultFile
  friend class TestUCDFaultFile; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  static
  void read(const char* filename,
	    const ALE::Obj<SieveMesh>& mesh,
	    const ALE::Obj<SieveMesh>& fault,
	    ALE::Obj<ALE::Mesh>& faultBd);

}; // UCDFaultFile

#endif // pylith_meshio_ucdfaultfile_hh


// End of file
