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

/** @file libsrc/faults/FaultCohesiveKin.hh
 *
 * @brief C++ abstract base class for a fault surface implemented with
 * cohesive elements.
 */

#if !defined(pylith_faults_faultcohesivekin_hh)
#define pylith_faults_faultcohesivekin_hh

#include "FaultCohesive.hh"

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class FaultCohesiveKin;
    class TestFaultCohesiveKin; // unit testing

    class EqKinSrc; // HOLDSA EqKinSrc
  } // faults
} // pylith

/// C++ implementation for a fault surface with kinematic (prescribed)
/// slip implemented with cohesive elements.
class pylith::faults::FaultCohesiveKin : public FaultCohesive
{ // class FaultCohesiveKin
  friend class TestFaultCohesiveKin; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesiveKin(void);

  /// Destructor.
  virtual
  ~FaultCohesiveKin(void);

  /** Create copy of fault.
   *
   * @returns Copy of fault.
   */
  Fault* clone(void) const;  

  /** Set kinematic earthquake source.
   *
   * @param src Kinematic earthquake source.
   */
  void eqsrc(EqKinSrc* src);

  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh PETSc mesh
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction).
   */
  void initialize(const ALE::Obj<ALE::Mesh>& mesh,
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
   * @param disp Displacement field
   * @param mesh Finite-element mesh
   */
  void setField(const ALE::Obj<real_section_type>& disp,
		const ALE::Obj<Mesh>& mesh);
  
  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param m Fault to copy
   */
  FaultCohesiveKin(const FaultCohesiveKin& m);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const FaultCohesiveKin& operator=(const FaultCohesiveKin& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  EqKinSrc* _eqsrc; ///< Kinematic earthquake source information

}; // class FaultCohesiveKin

#include "FaultCohesiveKin.icc" // inline methods

#endif // pylith_faults_faultcohesivekin_hh


// End of file 
