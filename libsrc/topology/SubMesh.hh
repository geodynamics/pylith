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
 * @file libsrc/topology/SubMesh.hh
 *
 * @brief C++ PyLith finite-element mesh.
 *
 * Extends Sieve mesh over subset of domain to include coordinate
 * system associated with domain.
 */

#if !defined(pylith_topology_submesh_hh)
#define pylith_topology_submesh_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations
#include "spatialdata/geocoords/geocoordsfwd.hh" // forward declarations

#define NEWPYLITHMESH 1 
#include "pylith/utils/sievetypes.hh"

// SubMesh -----------------------------------------------------------------
class pylith::topology::SubMesh
{ // SubMesh
  friend class TestSubMesh; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  SubMesh(void);

  /** Constructor with mesh and label for vertices marking boundary.
   *
   * @param mesh Finite-element mesh over domain.
   * @param label Label for vertices marking boundary.
   */
  SubMesh(const Mesh& mesh,
	  const char* label);

  /// Default destructor
  ~SubMesh(void);

  /** Create Sieve mesh.
   *
   * @param mesh Finite-element mesh over domain.
   * @param label Label for vertices marking boundary.
   */
  void createSubMesh(const Mesh& mesh,
		     const char* label); 

  /** Get Sieve mesh.
   *
   * @returns Sieve mesh.
   */
  const ALE::Obj<SieveSubMesh>& sieveMesh(void) const;

  /** Get Sieve mesh.
   *
   * @returns Sieve mesh.
   */
  ALE::Obj<SieveSubMesh>& sieveMesh(void);

  /** Get coordinate system.
   *
   * @returns Coordinate system.
   */
  const spatialdata::geocoords::CoordSys* coordsys(void) const;

  /** Set debug flag.
   *
   * @param value Turn on debugging if true.
   */
   void debug(const bool value);

  /** Get debug flag.
   *
   * @param Get debugging flag.
   */
   bool debug(void) const;

  /** Get dimension of mesh.
   *
   * @returns Dimension of mesh.
   */
  int dimension(void) const;

  /** Get MPI communicator associated with mesh.
   *
   * @returns MPI communicator.
   */
  const MPI_Comm comm(void) const;
    
  /// Initialize the finite-element mesh.
  void initialize(void);

  /** Print mesh to stdout.
   *
   * @param label Label for mesh.
   */
  void view(const char* label);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  ALE::Obj<SieveSubMesh> _mesh; ///< Sieve mesh.
  spatialdata::geocoords::CoordSys* _coordsys; ///< Coordinate system.
  bool _debug; ///< Debugging flag for mesh.
  
// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  SubMesh(const SubMesh&); ///< Not implemented
  const SubMesh& operator=(const SubMesh&); ///< Not implemented

}; // SubMesh

#include "SubMesh.icc"

#endif // pylith_topology_submesh_hh


// End of file
