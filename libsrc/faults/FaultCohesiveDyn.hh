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

#include "FaultCohesive.hh"

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class FaultCohesiveDyn;
    class TestFaultCohesiveDyn; // unit testing
  } // faults
} // pylith

/// @brief C++ implementation for a fault surface with spontaneous
/// (dynamic) slip implemented with cohesive elements.
class pylith::faults::FaultCohesiveDyn : public FaultCohesive
{ // class FaultCohesiveDyn
  friend class TestFaultCohesiveDyn; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesiveDyn(void);

  /// Destructor.
  virtual
  ~FaultCohesiveDyn(void);

  /** Create copy of fault.
   *
   * @returns Copy of fault.
   */
  Fault* clone(void) const;  

  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh PETSc mesh
   * @param cs Coordinate system for mesh
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction).
   */
  void initialize(const ALE::Obj<ALE::Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const double_array& upDir);

  /** Integrate contribution of cohesive cells to residual term.
   *
   * @param residual Residual field (output)
   * @param disp Displacement field at time t
   * @param mesh Finite-element mesh
   */
  void integrateResidual(const ALE::Obj<real_section_type>& residual,
			 const ALE::Obj<real_section_type>& disp,
			 const ALE::Obj<Mesh>& mesh);

  /** Compute Jacobian matrix (A) associated with operator.
   *
   * @param mat Sparse matrix
   * @param disp Displacement field
   * @param mesh Finite-element mesh
   */
  void integrateJacobian(PetscMat* mat,
			 const ALE::Obj<real_section_type>& dispT,
			 const ALE::Obj<Mesh>& mesh);
  
  /** Set field.
   *
   * @param t Current time
   * @param disp Displacement field
   * @param mesh Finite-element mesh
   */
  void setField(const double t,
		const ALE::Obj<real_section_type>& disp,
		const ALE::Obj<Mesh>& mesh);
  
  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param m Fault to copy
   */
  FaultCohesiveDyn(const FaultCohesiveDyn& m);

  /** Cohesive cells use Lagrange multiplier constraints?
   *
   * @returns True if implementation using Lagrange multiplier
   * constraints, false otherwise.
   */
  bool _useLagrangeConstraints(void) const;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

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
