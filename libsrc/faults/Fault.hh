// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/Fault.hh
 *
 * @brief C++ abstract base class for Fault object.
 *
 * Interface definition for fault.
 *
 * The fault id is associated with the material-id for the fault and
 * the label is associated with the group of vertices that define the
 * fault surface.
 */

#if !defined(pylith_faults_fault_hh)
#define pylith_faults_fault_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Field<SubMesh>, SubMesh
#include "pylith/utils/arrayfwd.hh" // USES double_array

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB

#include <string> // HASA std::string

// Fault ----------------------------------------------------------------
class pylith::faults::Fault
{ // class Fault
  friend class TestFault; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  Fault(void);

  /// Destructor.
  virtual
  ~Fault(void);

  /** Set material identifier of fault.
   *
   * @param value Fault identifier
   */
  void id(const int value);

  /** Get material identifier of fault.
   *
   * @returns Fault identifier
   */
  int id(void) const;

  /** Set label of group of vertices associated with fault.
   *
   * @param value Label of fault
   */
  void label(const char* value);

  /** Get label of group of vertices associated with fault.
   *
   * @returns Label of fault
   */
  const char* label(void) const;

  /** Adjust mesh topology for fault implementation.
   *
   * @param mesh PETSc mesh
   */
  virtual
  void adjustTopology(const topology::Mesh& mesh,
		      const bool flipFault =false) = 0;

  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh PETSc mesh
   * @param cs Coordinate system for mesh
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction; only applies to fault surfaces in a 3-D domain).
   * @param normalDir General preferred direction for fault normal
   *   (used to pick which of two possible normal directions for
   *   interface; only applies to fault surfaces in a 3-D domain).
   * @param matDB Database of bulk elastic properties for fault region
   *   (used to improve conditioning of Jacobian matrix)
   */
  virtual
  void initialize(const topology::Mesh& mesh,
		  const double_array& upDir,
		  const double_array& normalDir,
		  spatialdata::spatialdb::SpatialDB* matDB) = 0;

  /** Get mesh associated with fault fields.
   *
   * @returns PETSc mesh object
   */
  const topology::SubMesh& faultMesh(void) const;

  /** Get vertex field associated with integrator.
   *
   * @param name Name of vertex field.
   * @param fields Solution fields.
   * @returns Vertex field.
   */
  virtual
  const topology::Field<topology::SubMesh>&
  vertexField(const char* name,
	      const topology::SolutionFields& fields) = 0;

  /** Get cell field associated with integrator.
   *
   * @param name Name of cell field.
   * @param fields Solution fields.
   * @returns Cell field.
   */
  virtual
  const topology::Field<topology::SubMesh>&
  cellField(const char* name,
	    topology::SolutionFields& fields) = 0;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  Fault(const Fault&); ///< Not implemented
  const Fault& operator=(const Fault&); ///< Not implemented

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  topology::SubMesh* _faultMesh; ///< Mesh over fault surface

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  int _id; ///< Fault identifier
  std::string _label; ///< Label of fault

}; // class Fault

#include "Fault.icc" // inline methods

#endif // pylith_faults_fault_hh


// End of file 
