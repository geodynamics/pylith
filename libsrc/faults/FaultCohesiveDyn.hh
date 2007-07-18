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

/** @file libsrc/faults/FaultCohesiveDyn.hh
 *
 * @brief C++ implementation for a fault surface with spontaneous
 * (dynamic) slip implemented with cohesive elements.
 *
 * The ordering of vertices in a cohesive cell is the vertices on the
 * POSITIVE/NEGATIVE (CHECK WHICH IT IS) side of the fault and then the
 * corresponding entries on the other side of the fault.
 */

#if !defined(pylith_faults_faultcohesivedyn_hh)
#define pylith_faults_faultcohesivedyn_hh

#include "FaultCohesive.hh" // ISA FaultCohesive
#include "pylith/feassemble/Integrator.hh" // ISA Constraint

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class FaultCohesiveDyn;
    class TestFaultCohesiveDyn; // unit testing
  } // faults
} // pylith

/// @brief C++ implementation for a fault surface with spontaneous
/// (dynamic) slip implemented with cohesive elements.
class pylith::faults::FaultCohesiveDyn : public FaultCohesive,
					 public feassemble::Integrator
{ // class FaultCohesiveDyn
  friend class TestFaultCohesiveDyn; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesiveDyn(void);

  /// Destructor.
  virtual
  ~FaultCohesiveDyn(void);

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
  void initialize(const ALE::Obj<ALE::Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const double_array& upDir,
		  const double_array& normalDir,
		  spatialdata::spatialdb::SpatialDB* matDB);

  /** Integrate contribution of cohesive cells to residual term.
   *
   * @param residual Residual field (output)
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void integrateResidual(const ALE::Obj<real_section_type>& residual,
			 const double t,
			 topology::FieldsManager* const fields,
			 const ALE::Obj<Mesh>& mesh);

  /** Compute Jacobian matrix (A) associated with operator.
   *
   * @param mat Sparse matrix
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void integrateJacobian(PetscMat* mat,
			 const double t,
			 topology::FieldsManager* const fields,
			 const ALE::Obj<Mesh>& mesh);
  
  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const ALE::Obj<Mesh>& mesh);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Cohesive cells use Lagrange multiplier constraints?
   *
   * @returns True if implementation using Lagrange multiplier
   * constraints, false otherwise.
   */
  bool _useLagrangeConstraints(void) const;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesiveDyn(const FaultCohesiveDyn& m);

  /// Not implemented
  const FaultCohesiveDyn& operator=(const FaultCohesiveDyn& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Orientation of fault surface at vertices (fiber dimension is
  /// nonzero only at constraint vertices)
  ALE::Obj<real_section_type> _orientation;

}; // class FaultCohesiveDyn

#include "FaultCohesiveDyn.icc" // inline methods

#endif // pylith_faults_faultcohesivedyn_hh


// End of file 
