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
#include "meshiofwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Mesh, SubMesh

#include "pylith/utils/sievetypes.hh" // USES ALE::Obj, ALE::Mesh

// UCDFaultFile ---------------------------------------------------------
class pylith::meshio::UCDFaultFile
{ // UCDFaultFile
  friend class TestUCDFaultFile; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  static
  void read(const char* filename,
	    topology::SubMesh* faultMesh,
	    ALE::Obj<ALE::Mesh>& faultBoundary,
	    const topology::Mesh& mesh);

}; // UCDFaultFile

#endif // pylith_meshio_ucdfaultfile_hh


// End of file
