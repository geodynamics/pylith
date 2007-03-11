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
 * @file pylith/feassemble/IntegratorExplicit.hh
 *
 * @brief Abstract base class for explicit time integration of
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

#if !defined(pylith_feassemble_integratorexplicit_hh)
#define pylith_feassemble_integratorexplicit_hh

#include "pylith/utils/petscfwd.h" // USES PetscMat

#include "Integrator.hh" // ISA Integrator

namespace pylith {
  namespace feassemble {
    class IntegratorExplicit;
    class TestIntegratorExplicit;
  } // feassemble
} // pylith

class pylith::feassemble::IntegratorExplicit : public Integrator
{ // Integrator
  friend class TestIntegratorExplicit; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  IntegratorExplicit(void);

  /// Destructor
  virtual
  ~IntegratorExplicit(void);

  /// Create a copy of this object.
  virtual
  IntegratorExplicit* clone(void) const = 0;

  /** Set quadrature for integrating finite-element quantities.
   *
   * @param q Quadrature for integrating.
   */
  void quadrature(const Quadrature* q);

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

  /** Integrate residual term (b) for dynamic elasticity term 
   * for finite elements.
   *
   * @param residual Residual field (output)
   * @param fieldInT Input field at time t
   * @param fieldInTmdt Input field at time t-dt
   * @param coordinates Field of cell vertex coordinates
   */
  virtual 
  void integrateResidual(const ALE::Obj<real_section_type>& residual,
			 const ALE::Obj<real_section_type>& fieldInT,
			 const ALE::Obj<real_section_type>& fieldInTmdt,
			 const ALE::Obj<real_section_type>& coordinates) = 0;

  /** Compute matrix (A) associated with operator.
   *
   * @param mat Sparse matrix
   * @param fieldIn Input field at time t
   * @param coordinates Field of cell vertex coordinates
   */
  virtual 
  void integrateJacobian(PetscMat* mat,
			 const ALE::Obj<real_section_type>& fieldIn,
			 const ALE::Obj<real_section_type>& coordinates) = 0;
  
  /** Compute field (A) associated with lumped operator.
   *
   * @param fieldOut Output Jacobian field
   * @param fieldInT Input field at time t
   * @param coordinates Field of cell vertex coordinates
   */
  virtual 
  void integrateJacobian(const ALE::Obj<real_section_type>& fieldOut,
			 const ALE::Obj<real_section_type>& fieldInT,
			 const ALE::Obj<real_section_type>& coordinates) = 0;
  
// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param i Integrator to copy
   */
  IntegratorExplicit(const IntegratorExplicit& i);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const IntegratorExplicit& operator=(const IntegratorExplicit&);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  double _dt; ///< Time step for t -> t+dt
  double _dtm1; ///< Time step for t-dt1 -> t

}; // IntegratorExplicit

#endif // pylith_feassemble_integratorexplicit_hh


// End of file 
