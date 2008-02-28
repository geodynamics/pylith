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

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh, real_section_type
#include "pylith/utils/arrayfwd.hh" // USES double_array
#include "pylith/utils/vectorfields.hh" // USES VectorFieldEnum

#include <string> // HASA std::string

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class Fault;
    class TestFault; // unit testing
  } // faults

  namespace topology {
    class FieldsManager;
  } // topology
} // pylith

/// Namespace for spatialdata package
namespace spatialdata {
  namespace geocoords {
    class CoordSys;
  } // geocoords

  namespace spatialdb {
    class SpatialDB; // USES SpatialDB
  } // spatialdb
} // spatialdata

/// C++ abstract base class for Fault object.
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
  const std::string& label(void) const;

  /** Adjust mesh topology for fault implementation.
   *
   * @param mesh PETSc mesh
   */
  virtual
  void adjustTopology(const ALE::Obj<ALE::Mesh>& mesh) = 0;

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
  void initialize(const ALE::Obj<ALE::Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const double_array& upDir,
		  const double_array& normalDir,
		  spatialdata::spatialdb::SpatialDB* matDB) = 0;

  /** Get mesh associated with fault fields.
   *
   * @returns PETSc mesh object
   */
  const ALE::Obj<ALE::Mesh>& faultMesh(void) const;

  /** Get vertex field associated with integrator.
   *
   * @param fieldType Type of field.
   * @param name Name of vertex field.
   * @param mesh PETSc mesh for problem.
   * @param fields Fields manager.
   * @returns Vertex field.
   */
  virtual
  const ALE::Obj<real_section_type>&
  vertexField(VectorFieldEnum* fieldType,
	      const char* name,
	      const ALE::Obj<Mesh>& mesh,
	      topology::FieldsManager* const fields) = 0;

  /** Get cell field associated with integrator.
   *
   * @param fieldType Type of field.
   * @param name Name of cell field.
   * @param mesh PETSc mesh for problem.
   * @param fields Fields manager.
   * @returns Cell field.
   */
  virtual
  const ALE::Obj<real_section_type>&
  cellField(VectorFieldEnum* fieldType,
	    const char* name,
	    const ALE::Obj<Mesh>& mesh,
	    topology::FieldsManager* const fields) = 0;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  Fault(const Fault& m);

  /// Not implemented
  const Fault& operator=(const Fault& m);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  ALE::Obj<ALE::Mesh> _faultMesh; ///< Mesh over fault surface

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  int _id; ///< Fault identifier
  std::string _label; ///< Label of fault

}; // class Fault

#include "Fault.icc" // inline methods

#endif // pylith_faults_fault_hh


// End of file 
