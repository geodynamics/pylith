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
 * @file pylith/feassemble/IntegratorImplicit.hh
 *
 * @brief Abstract base class for implicit time integration of
 * finite-element actions.
 *
 * Note: This object operates on a single finite-element family, which
 * is defined by the quadrature and a database of material property
 * parameters.
 *
 * Computes terms A and b in A(t) u(t+dt) = b(u(t), u(t-dt)), where
 * A(t) is a sparse matrix or vector, u(t+dt) is the field we want to
 * compute at time t+dt and b is a vector that depends on the field at
 * time t and t-dt.
 */

#if !defined(pylith_feassemble_integratorimplicit_hh)
#define pylith_feassemble_integratorimplicit_hh

#include "pylith/utils/petscfwd.h" // USES PetscMat

#include "Integrator.hh" // ISA Integrator

namespace pylith {
  namespace feassemble {
    class IntegratorImplicit;
    class TestIntegratorImplicit;
  } // feassemble
} // pylith

class pylith::feassemble::IntegratorImplicit : public Integrator
{ // Integrator
  friend class TestIntegratorImplicit; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  IntegratorImplicit(void);

  /// Destructor
  virtual
  ~IntegratorImplicit(void);

  /// Create a copy of this object.
  virtual
  IntegratorImplicit* clone(void) const = 0;

  /** Set time step for advancing from time t to time t+dt.
   *
   * @param dt Time step
   */
  void timeStep(const double dt);

  /** Get stable time step for advancing from time t to time t+dt.
   *
   * Default is current time step.
   *
   * @returns Time step
   */
  virtual
  double stableTimeStep(void) const;

  /** Integrate constant term (b) for dynamic elasticity term 
   * for finite elements.
   *
   * @param fieldOut Constant field (output)
   * @param fieldInT Input field at time t
   * @param mesh Mesh object
   */
  virtual 
  void integrateConstant(const ALE::Obj<real_section_type>& fieldOut,
			 const ALE::Obj<real_section_type>& fieldInT,
			 const ALE::Obj<Mesh>& mesh) = 0;

  /** Integrate residual for quasi-static finite elements.
   *
   * @param fieldOut Constant field (output)
   * @param fieldInT Input field at time t
   * @param mesh Mesh object
   */
  virtual 
  void integrateResidual(const ALE::Obj<real_section_type>& fieldOut,
			 const ALE::Obj<real_section_type>& fieldInT,
			 const ALE::Obj<Mesh>& mesh) = 0;

  /** Compute Jacobian matrix (A) associated with operator.
   *
   * @param mat Sparse matrix
   * @param fieldIn Input field at time t
   * @param coordinates Field of cell vertex coordinates
   */
  virtual 
  void integrateJacobian(PetscMat* mat,
			 const ALE::Obj<real_section_type>& fieldIn,
			 const ALE::Obj<Mesh>& mesh) = 0;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param i Integrator to copy
   */
  IntegratorImplicit(const IntegratorImplicit& i);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const IntegratorImplicit& operator=(const IntegratorImplicit&);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  double _dt; ///< Time step for t -> t+dt
  double _dtm1; ///< Time step for t-dt1 -> t

}; // IntegratorImplicit

#endif // pylith_feassemble_integratorimplicit_hh


// End of file 
