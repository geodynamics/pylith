// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/IntegratorElasticityLgDeform.hh
 *
 * @brief Object containing general elasticity operations for implicit
 * and explicit time integration of the elasticity equation for large
 * rigid body rotations and small strains.
 */

#if !defined(pylith_feassemble_integratorelasticitylgdeform_hh)
#define pylith_feassemble_integratorelasticitylgdeform_hh

// Include directives ---------------------------------------------------
#include "IntegratorElasticity.hh" // ISA IntegratorElasticity

// IntegratorElasticity -------------------------------------------------
/** @brief General elasticity operations for implicit and explicit
 * time integration of the elasticity equation for large rigid body
 * rotations and small strains.
 */
class pylith::feassemble::IntegratorElasticityLgDeform : public IntegratorElasticity
{ // IntegratorElasticityLgDeform
  friend class TestIntegratorElasticityLgDeform; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :

  typedef void (*totalStrain_fn_type)(scalar_array*,
				      const scalar_array&,
				      const int);
  

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  IntegratorElasticityLgDeform(void);

  /// Destructor
  virtual
  ~IntegratorElasticityLgDeform(void);

  /** Determine whether we need to recompute the Jacobian.
   *
   * @returns True if Jacobian needs to be recomputed, false otherwise.
   */
  bool needNewJacobian(void);

  /** Update state variables as needed.
   *
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void updateStateVars(const PylithScalar t,
		       topology::SolutionFields* const fields);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Calculate stress or strain field from solution field.
   *
   * @param field Field in which to store stress or strain.
   * @param name Name of field to compute ('total-strain' or 'stress')
   * @param fields Manager for solution fields.
   */
  void _calcStrainStressField(topology::Field* field,
			      const char* name,
			      topology::SolutionFields* const fields);

  /** Integrate elasticity term in residual for 2-D cells.
   *
   * @param stress Stress tensor for cell at quadrature points.
   * @param disp Displacement field at cell's DOF.
   */
  void _elasticityResidual2D(const scalar_array& stress,
			     const scalar_array& disp);

  /** Integrate elasticity term in residual for 3-D cells.
   *
   * @param stress Stress tensor for cell at quadrature points.
   * @param disp Displacement field at cell's DOF.
   */
  void _elasticityResidual3D(const scalar_array& stress,
			     const scalar_array& disp);

  /** Integrate elasticity term in Jacobian for 2-D cells.
   *
   * @param elasticConsts Matrix of elasticity constants at quadrature points.
   * @param stress Stress tensor for cell at quadrature points.
   * @param disp Displacement field at cell's DOF.
   */
  void _elasticityJacobian2D(const scalar_array& elasticConsts,
			     const scalar_array& stress,
			     const scalar_array& disp);

  /** Integrate elasticity term in Jacobian for 3-D cells.
   *
   * @param elasticConsts Matrix of elasticity constants at quadrature points.
   * @param stress Stress tensor for cell at quadrature points.
   * @param disp Displacement field at cell's DOF.
   */
  void _elasticityJacobian3D(const scalar_array& elasticConsts,
			     const scalar_array& stress,
			     const scalar_array& disp);

  /** Calculate Green-Lagrange strain tensor at quadrature points of a
   *  1-D cell.
   *
   * @param strain Green-Lagrange strain tensor at quadrature points (output).
   * @param deform Deformation tensor at quadrature points.
   * @param numQuadPts Number of quadrature points.
   */
  static
  void _calcTotalStrain1D(scalar_array* strain,
			  const scalar_array& deform,
			  const int numQuadPts);

  /** Calculate Green-Lagrange strain tensor at quadrature points of a
   *  2-D cell.
   *
   * @param strain Green-Lagrange strain tensor at quadrature points (output).
   * @param deform Deformation tensor at quadrature points.
   * @param numQuadPts Number of quadrature points.
   */
  static
  void _calcTotalStrain2D(scalar_array* strain,
			  const scalar_array& deform,
			  const int numQuadPts);

  /** Calculate Green-Lagrange strain tensor at quadrature points of a
   *  3-D cell.
   *
   * @param strain Green-Lagrange strain tensor at quadrature points (output).
   * @param deform Deformation tensor at quadrature points.
   * @param numQuadPts Number of quadrature points.
   */
  static
  void _calcTotalStrain3D(scalar_array* strain,
			  const scalar_array& deform,
			  const int numQuadPts);

  /** Calculate deformation tensor.
   *
   * @param deform Deformation tensor for cell at quadrature points.
   * @param basisDeriv Derivatives of basis functions at quadrature points.
   * @param vertices Coordinates of vertices in reference cell.
   * @param disp Displacements of DOF of cell.
   * @param numBasis Number of basis functions for cell.
   * @param numQuadPts Number of quadrature points.
   * @param dim Dimension of cell.
   */
  static
  void _calcDeformation(scalar_array* deform,
			const scalar_array& basisDeriv,
			const PylithScalar* vertices,
			const PylithScalar* disp,
			const int numBasis,
			const int numQuadPts,
			const int dim);
  
// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented.
  IntegratorElasticityLgDeform(const IntegratorElasticityLgDeform&);

  /// Not implemented
  const IntegratorElasticityLgDeform& operator=(const IntegratorElasticityLgDeform&);

}; // IntegratorElasticityLgDeform

#endif // pylith_feassemble_integratorelasticitylgdeform_hh


// End of file 
