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
 * @file pylith/feassemble/Integrator.hh
 *
 * @brief Abstract base class for integration of finite-element
 * actions.
 *
 * Note: Each object operates on a single finite-element family, which
 * is defined by the quadrature and a database of material property
 * parameters.
 */

#if !defined(pylith_feassemble_integrator_hh)
#define pylith_feassemble_integrator_hh

#include "pylith/utils/sievetypes.hh" // USES real_section_type
#include "pylith/utils/petscfwd.h" // USES PetscMat

namespace pylith {
  namespace feassemble {
    class Integrator;
    class TestIntegrator;

    class Quadrature; // HOLDSA Quadrature
  } // feassemble

  namespace topology {
    class FieldsManager;
  } // topology
} // pylith

class pylith::feassemble::Integrator
{ // Integrator
  friend class TestIntegrator; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Integrator(void);

  /// Destructor
  virtual
  ~Integrator(void);

  /** Set quadrature for integrating finite-element
   * quantities. Quadrature should already be initialized.
   *
   * @param q Quadrature for integrating.
   */
  void quadrature(const Quadrature* q);

  /** Set time step for advancing from time t to time t+dt.
   *
   * @param dt Time step
   */
  virtual
  void timeStep(const double dt);

  /** Get stable time step for advancing from time t to time t+dt.
   *
   * Default is current time step.
   *
   * @returns Time step
   */
  virtual
  double stableTimeStep(void) const;

  /** Check whether Jacobian needs to be recomputed.
   *
   * @returns True if Jacobian needs to be recomputed, false otherwise.
   */
  virtual
  bool needNewJacobian(void) const;

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  virtual 
  void integrateResidual(const ALE::Obj<real_section_type>& residual,
			 topology::FieldsManager* const fields,
			 const ALE::Obj<Mesh>& mesh) = 0;

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param mat Sparse matrix
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  virtual 
  void integrateJacobian(PetscMat* mat,
			 topology::FieldsManager* const fields,
			 const ALE::Obj<Mesh>& mesh) = 0;

  /** Update state variables as needed.
   *
   * @param field Current solution field.
   * @param mesh Finite-element mesh
   */
  virtual
  void updateState(const ALE::Obj<real_section_type>& field,
		   const ALE::Obj<Mesh>& mesh);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /// Initialize vector containing result of integration action for cell.
  void _initCellVector(void);

  /// Zero out vector containing result of integration actions for cell.
  void _resetCellVector(void);

  /// Initialize matrix containing result of integration for cell.
  void _initCellMatrix(void);

  /// Zero out matrix containing result of integration for cell.
  void _resetCellMatrix(void);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  // Not implemented.
  Integrator(const Integrator& i);

  /// Not implemented
  const Integrator& operator=(const Integrator&);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  double _dt; ///< Time step for t -> t+dt

  Quadrature* _quadrature; ///< Quadrature for integrating finite-element

  /// Vector local to cell containing result of integration action
  real_section_type::value_type* _cellVector;

  /// Matrix local to cell containing result of integration
  real_section_type::value_type* _cellMatrix;

  /// True if we need to recompute Jacobian for operator, false otherwise.
  /// Default is false;
  bool _needNewJacobian;

}; // Integrator

#endif // pylith_feassemble_integrator_hh

// End of file 
