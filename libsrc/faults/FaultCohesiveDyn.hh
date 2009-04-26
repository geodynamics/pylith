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

// Include directives ---------------------------------------------------
#include "FaultCohesive.hh" // ISA FaultCohesive

#include "pylith/topology/SubMesh.hh" // ISA Integrator<Quadrature<SubMesh> >
#include "pylith/feassemble/Quadrature.hh" // ISA Integrator<Quadrature>
#include "pylith/feassemble/Integrator.hh" // ISA Integrator

// FaultCohesiveDyn -----------------------------------------------------
class pylith::faults::FaultCohesiveDyn : public FaultCohesive,
					 public feassemble::Integrator<feassemble::Quadrature<topology::SubMesh> >
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
   * @param mesh Finite-element mesh.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction; only applies to fault surfaces in a 3-D domain).
   * @param normalDir General preferred direction for fault normal
   *   (used to pick which of two possible normal directions for
   *   interface; only applies to fault surfaces in a 3-D domain).
   * @param matDB Database of bulk elastic properties for fault region
   *   (used to improve conditioning of Jacobian matrix)
   */
  void initialize(const topology::Mesh& mesh,
		  const double upDir[3],
		  const double normalDir[3],
		  spatialdata::spatialdb::SpatialDB* matDB);

  /** Integrate contribution of cohesive cells to residual term.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidual(const topology::Field<topology::Mesh>& residual,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateJacobian(topology::Jacobian* jacobian,
			 const double t,
			 topology::SolutionFields* const fields);
  
  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const topology::Mesh& mesh) const;

  /** Get vertex field associated with integrator.
   *
   * @param name Name of vertex field.
   * @param fields Solution fields.
   *
   * @returns Vertex field.
   */
  const topology::Field<topology::SubMesh>&
  vertexField(const char* name,
	      const topology::SolutionFields& fields);
  
  /** Get cell field associated with integrator.
   *
   * @param name Name of cell field.
   * @param fields Solution fields.
   *
   * @returns Cell field.
   */
  const topology::Field<topology::SubMesh>&
  cellField(const char* name,
	    const topology::SolutionFields& fields);

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
  FaultCohesiveDyn(const FaultCohesiveDyn&);

  /// Not implemented
  const FaultCohesiveDyn& operator=(const FaultCohesiveDyn&);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Orientation of fault surface at vertices (fiber dimension is
  /// nonzero only at constraint vertices)
  topology::Field<topology::SubMesh>*  _orientation;

}; // class FaultCohesiveDyn

#include "FaultCohesiveDyn.icc" // inline methods

#endif // pylith_faults_faultcohesivedyn_hh


// End of file 
