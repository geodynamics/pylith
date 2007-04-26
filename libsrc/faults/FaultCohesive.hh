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

/** @file libsrc/faults/FaultCohesive.hh
 *
 * @brief C++ abstract base class for a fault surface implemented with
 * cohesive elements.
 */

#if !defined(pylith_faults_faultcohesive_hh)
#define pylith_faults_faultcohesive_hh

#include "Fault.hh" // ISA Fault
#include "pylith/feassemble/Integrator.hh" // ISA Integrator

#include "pylith/utils/sievefwd.hh" // HOLDSA PETSc Mesh
#include "pylith/utils/petscfwd.h" // USES PetscMat

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class FaultCohesive;
    class TestFaultCohesive; // unit testing
  } // faults
} // pylith

/// C++ abstract base class for a fault surface implemented with
/// cohesive elements.
class pylith::faults::FaultCohesive : public Fault, 
				      public feassemble::Integrator
{ // class FaultCohesive
  friend class TestFaultCohesive; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesive(void);

  /// Destructor.
  virtual
  ~FaultCohesive(void);

  /** Adjust mesh topology for fault implementation.
   *
   * @param mesh PETSc mesh
   */
  void adjustTopology(const ALE::Obj<ALE::Mesh>& mesh);

  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh PETSc mesh
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction).
   */
  virtual
  void initialize(const ALE::Obj<ALE::Mesh>& mesh,
		  const double_array& upDir) = 0;

  /** Integrate contribution of cohesive cells to residual term.
   *
   * @param residual Residual field (output)
   * @param disp Displacement field at time t
   * @param mesh Finite-element mesh
   */
  virtual
  void integrateResidual(const ALE::Obj<real_section_type>& residual,
			 const ALE::Obj<real_section_type>& disp,
			 const ALE::Obj<Mesh>& mesh) = 0;

  /** Compute Jacobian matrix (A) associated with operator.
   *
   * @param mat Sparse matrix
   * @param disp Displacement field
   * @param mesh Finite-element mesh
   */
  virtual
  void integrateJacobian(PetscMat* mat,
			 const ALE::Obj<real_section_type>& dispT,
			 const ALE::Obj<Mesh>& mesh) = 0;
  
  /** Set field.
   *
   * @param disp Displacement field
   * @param mesh Finite-element mesh
   */
  virtual
  void setField(const ALE::Obj<real_section_type>& disp,
		const ALE::Obj<Mesh>& mesh) = 0;
  
  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param m Fault to copy
   */
  FaultCohesive(const FaultCohesive& m);

  /** Compute weighted orientation of fault for cohesive cell between
   * 1-D elements. Orientation is either at vertices or quadrature
   * points, depending on whether the arguments have been evaluated at
   * the vertices or quadrature points.
   *
   * The orientation is returned as an array of direction cosines.
   *
   * size = numPts*spaceDim*spaceDim
   * index = iPt*spaceDim*spaceDim + iDir*spaceDim + iComponent
   *
   * @param orientation Array of direction cosines.
   * @param jacobian Jacobian matrix at point.
   * @param jacobianDet Determinant of Jacobian matrix at point.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction).
   */
  void orient1D(double_array* orientation,
		const double_array& jacobian,
		const double_array& jacobianDet,
		const double_array& upDir);
		
  /** Compute weighted orientation of fault for cohesive cell between
   * 2-D elements. Orientation is either at vertices or quadrature
   * points, depending on whether the arguments have been evaluated at
   * the vertices or quadrature points.
   *
   * The orientation is returned as an array of direction cosines.
   *
   * size = numPts*spaceDim*spaceDim
   * index = iPt*spaceDim*spaceDim + iDir*spaceDim + iComponent
   *
   * @param orientation Array of direction cosines.
   * @param jacobian Jacobian matrix at point.
   * @param jacobianDet Determinant of Jacobian matrix at point.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction).
   */
  void orient2D(double_array* orientation,
		const double_array& jacobian,
		const double_array& jacobianDet,
		const double_array& upDir);
		
  /** Compute weighted orientation of fault for cohesive cell between
   * 3-D elements. Orientation is either at vertices or quadrature
   * points, depending on whether the arguments have been evaluated at
   * the vertices or quadrature points.
   *
   * The orientation is returned as an array of direction cosines.
   *
   * size = numPts*spaceDim*spaceDim
   * index = iPt*spaceDim*spaceDim + iDir*spaceDim + iComponent
   *
   * @param orientation Array of direction cosines.
   * @param jacobian Jacobian matrix at point.
   * @param jacobianDet Determinant of Jacobian matrix at point.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction).
   */
  void orient3D(double_array* orientation,
		const double_array& jacobian,
		const double_array& jacobianDet,
		const double_array& upDir);
		

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const FaultCohesive& operator=(const FaultCohesive& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  ALE::Obj<ALE::Mesh>* _faultMesh;

}; // class FaultCohesive

#endif // pylith_faults_faultcohesive_hh


// End of file 
