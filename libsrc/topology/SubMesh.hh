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

// SubMesh -----------------------------------------------------------------
template<typename mesh_type>
class pylith::topology::SubMesh
{ // SubMesh
  friend class TestSubMesh; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public:

  // Typedefs for basic types associated with Sieve mesh.
  // SieveMesh, RealSection, and IntSection are used in templated code.
  typedef typename mesh_type::SieveSubMesh SieveMesh;
  typedef typename mesh_type::RealSection  RealSection;
  typedef typename mesh_type::IntSection IntSection;

  // Sieve mesh for higher level domain (mesh, not submesh)
  typedef typename mesh_type::SieveMesh DomainSieveMesh;

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  SubMesh(void);

  /** Constructor with mesh and label for vertices marking boundary.
   *
   * @param mesh Finite-element mesh over domain.
   * @param label Label for vertices marking boundary.
   */
  SubMesh(const mesh_type& mesh,
	  const char* label);

  /// Default destructor
  ~SubMesh(void);

  /** Create Sieve mesh.
   *
   * @param mesh Finite-element mesh over domain.
   * @param label Label for vertices marking boundary.
   */
  void createSubMesh(const mesh_type& mesh,
		     const char* label); 

  /** Get Sieve mesh.
   *
   * @returns Sieve mesh.
   */
  const ALE::Obj<SieveMesh>& sieveMesh(void) const;

  /** Get Sieve mesh.
   *
   * @returns Sieve mesh.
   */
  ALE::Obj<SieveMesh>& sieveMesh(void);

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

  ALE::Obj<SieveMesh> _mesh; ///< Sieve mesh.
  spatialdata::geocoords::CoordSys* _coordsys; ///< Coordinate system.
  bool _debug; ///< Debugging flag for mesh.
  
// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  SubMesh(const SubMesh&); ///< Not implemented
  const SubMesh& operator=(const SubMesh&); ///< Not implemented

}; // SubMesh

#include "SubMesh.icc"
#include "SubMesh.cc"

#endif // pylith_topology_submesh_hh


// End of file
