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
 * @file libsrc/feassemble/IntegratorElasticity.hh
 *
 * @brief Object containing general elasticity operations for implicit
 * and explicit time integration of the elasticity equation.
 */

#if !defined(pylith_feassemble_integratorelasticity_hh)
#define pylith_feassemble_integratorelasticity_hh

// Include directives ---------------------------------------------------
#include "feassemblefwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // HOLDSA Field
#include "pylith/materials/materialsfwd.hh" // HOLDSA Material

#include "Integrator.hh" // ISA Integrator

#include "pylith/utils/arrayfwd.hh" // USES std::vector, scalar_array

// IntegratorElasticity -------------------------------------------------
/** @brief General elasticity operations for implicit and explicit
 * time integration of the elasticity equation.
 */
class pylith::feassemble::IntegratorElasticity : public Integrator
{ // IntegratorElasticity
  friend class TestIntegratorElasticity; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :

  typedef void (*totalStrain_fn_type)(scalar_array*,
				      const scalar_array&,
				      const PylithScalar*,
				      const int,
				      const int,
				      const int);
  

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  IntegratorElasticity(void);

  /// Destructor
  virtual
  ~IntegratorElasticity(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set material.
   *
   * @param m Elastic material.
   */
  void material(materials::ElasticMaterial* m);

  /** Determine whether we need to recompute the Jacobian.
   *
   * @returns True if Jacobian needs to be recomputed, false otherwise.
   */
  virtual
  bool needNewJacobian(void);

  /** Initialize integrator.
   *
   * @param mesh Finite-element mesh.
   */
  void initialize(const topology::Mesh& mesh);
  
  /** Update state variables as needed.
   *
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  virtual
  void updateStateVars(const PylithScalar t,
		       topology::SolutionFields* const fields);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  virtual
  void verifyConfiguration(const topology::Mesh& mesh) const;

  /** Get output fields.
   *
   * @returns Output (buffer) fields.
   */
  const topology::Fields* outputFields(void) const;

  /** Get cell field associated with integrator.
   *
   * @param name Name of vertex field.
   * @param mesh Finite-element mesh for problem.
   * @param fields Fields manager.
   * @returns Cell field.
   */
  const topology::Field& cellField(const char* name,
				   const topology::Mesh& mesh,
				   topology::SolutionFields* const fields =0);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /// Initialize logger.
  void _initializeLogger(void);

  /** Allocate buffer for tensor field at quadrature points.
   *
   * @param mesh Finite-element mesh.
   */
  void _allocateTensorField(const topology::Mesh& mesh);

  /** Calculate stress or strain field from solution field.
   *
   * @param field Field in which to store stress or strain.
   * @param name Name of field to compute ('total-strain' or 'stress')
   * @param fields Manager for solution fields.
   */
  virtual
  void _calcStrainStressField(topology::Field* field,
			      const char* name,
			      topology::SolutionFields* const fields);

  /** Calculate stress field from total strain field. Stress field
   * replaces strain field in section.
   *
   * @param field Field in which to store stress.
   */
  virtual
  void _calcStressFromStrain(topology::Field* field);

  /** Integrate elasticity term in residual for 2-D cells.
   *
   * @param stress Stress tensor for cell at quadrature points.
   */
  virtual
  void _elasticityResidual2D(const scalar_array& stress);

  /** Integrate elasticity term in residual for 3-D cells.
   *
   * @param stress Stress tensor for cell at quadrature points.
   */
  virtual
  void _elasticityResidual3D(const scalar_array& stress);

  /** Integrate elasticity term in Jacobian for 2-D cells.
   *
   * @param elasticConsts Matrix of elasticity constants at quadrature points.
   */
  virtual
  void _elasticityJacobian2D(const scalar_array& elasticConsts);

  /** Integrate elasticity term in Jacobian for 3-D cells.
   *
   * @param elasticConsts Matrix of elasticity constants at quadrature points.
   */
  virtual
  void _elasticityJacobian3D(const scalar_array& elasticConsts);

  /** Compute total strain in at quadrature points of a cell.
   *
   * @param strain Strain tensor at quadrature points.
   * @param basisDeriv Derivatives of basis functions at quadrature points.
   * @param disp Displacement at vertices of cell.
   * @param numBasis Number of basis functions for cell.
   * @param spaceDim Spatial dimension.
   * @param numQuadPts Number of quadrature points.
   */
  static
  void _calcTotalStrain2D(scalar_array* strain,
			  const scalar_array& basisDeriv,
			  const PylithScalar* disp,
			  const int numBasis,
			  const int spaceDim,
			  const int numQuadPts);

  /** Compute total strain in at quadrature points of a cell.
   *
   * @param strain Strain tensor at quadrature points.
   * @param basisDeriv Derivatives of basis functions at quadrature points.
   * @param disp Displacement at vertices of cell.
   * @param numBasis Number of basis functions for cell.
   * @param spaceDim Spatial dimension.
   * @param numQuadPts Number of quadrature points.
   */
  static
  void _calcTotalStrain3D(scalar_array* strain,
			  const scalar_array& basisDeriv,
			  const PylithScalar* disp,
			  const int numBasis,
			  const int spaceDim,
			  const int numQuadPts);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  
  materials::ElasticMaterial* _material; ///< Material associated with integrator.

  topology::StratumIS* _materialIS; ///< Index set for material cells.
  
  topology::Fields* _outputFields; ///< Buffers for output.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented.
  IntegratorElasticity(const IntegratorElasticity&);

  /// Not implemented
  const IntegratorElasticity& operator=(const IntegratorElasticity&);

}; // IntegratorElasticity

#endif // pylith_feassemble_integratorelasticity_hh


// End of file 
