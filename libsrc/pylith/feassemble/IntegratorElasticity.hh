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
// Copyright (c) 2010-2012 University of California, Davis
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

#include "pylith/topology/Mesh.hh" // ISA Integrator<Mesh>
#include "Integrator.hh" // ISA Integrator

#include "pylith/utils/arrayfwd.hh" // USES std::vector, double_array

// IntegratorElasticity -------------------------------------------------
/** @brief General elasticity operations for implicit and explicit
 * time integration of the elasticity equation.
 */
class pylith::feassemble::IntegratorElasticity :
  public Integrator<Quadrature<topology::Mesh> >
{ // IntegratorElasticity
  friend class TestIntegratorElasticity; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :

  typedef void (*totalStrain_fn_type)(double_array*,
				      const double_array&,
				      const double_array&,
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
  bool needNewJacobian(void);

  /** Set flag for setting constraints for total field solution or
   *  incremental field solution.
   *
   * @param flag True if using incremental solution, false otherwise.
   */
  void useSolnIncr(const bool flag);

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
  void updateStateVars(const double t,
		       topology::SolutionFields* const fields);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const topology::Mesh& mesh) const;

  /** Get output fields.
   *
   * @returns Output (buffer) fields.
   */
  const topology::Fields<topology::Field<topology::Mesh> >*
  outputFields(void) const;

  /** Get cell field associated with integrator.
   *
   * @param name Name of vertex field.
   * @param mesh Finite-element mesh for problem.
   * @param fields Fields manager.
   * @returns Cell field.
   */
  const topology::Field<topology::Mesh>&
  cellField(const char* name,
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
  void _calcStrainStressField(topology::Field<topology::Mesh>* field,
			      const char* name,
			      topology::SolutionFields* const fields);

  /** Calculate stress field from total strain field. Stress field
   * replaces strain field in section.
   *
   * @param field Field in which to store stress.
   */
  virtual
  void _calcStressFromStrain(topology::Field<topology::Mesh>* field);
			      

  /** Integrate elasticity term in residual for 1-D cells.
   *
   * @param stress Stress tensor for cell at quadrature points.
   */
  void _elasticityResidual1D(const double_array& stress);

  /** Integrate elasticity term in residual for 2-D cells.
   *
   * @param stress Stress tensor for cell at quadrature points.
   */
  void _elasticityResidual2D(const double_array& stress);

  /** Integrate elasticity term in residual for 3-D cells.
   *
   * @param stress Stress tensor for cell at quadrature points.
   */
  void _elasticityResidual3D(const double_array& stress);

  /** Integrate elasticity term in Jacobian for 1-D cells.
   *
   * @param elasticConsts Matrix of elasticity constants at quadrature points.
   */
  void _elasticityJacobian1D(const double_array& elasticConsts);

  /** Integrate elasticity term in Jacobian for 2-D cells.
   *
   * @param elasticConsts Matrix of elasticity constants at quadrature points.
   */
  void _elasticityJacobian2D(const double_array& elasticConsts);

  /** Integrate elasticity term in Jacobian for 3-D cells.
   *
   * @param elasticConsts Matrix of elasticity constants at quadrature points.
   */
  void _elasticityJacobian3D(const double_array& elasticConsts);

  /** Compute total strain in at quadrature points of a cell.
   *
   * @param strain Strain tensor at quadrature points.
   * @param basisDeriv Derivatives of basis functions at quadrature points.
   * @param disp Displacement at vertices of cell.
   * @param dimension Dimension of cell.
   * @param numBasis Number of basis functions for cell.
   * @param numQuadPts Number of quadrature points.
   */
  static
  void _calcTotalStrain1D(double_array* strain,
			  const double_array& basisDeriv,
			  const double_array& disp,
			  const int numBasis,
			  const int numQuadPts);

  /** Compute total strain in at quadrature points of a cell.
   *
   * @param strain Strain tensor at quadrature points.
   * @param basisDeriv Derivatives of basis functions at quadrature points.
   * @param disp Displacement at vertices of cell.
   * @param numBasis Number of basis functions for cell.
   * @param numQuadPts Number of quadrature points.
   */
  static
  void _calcTotalStrain2D(double_array* strain,
			  const double_array& basisDeriv,
			  const double_array& disp,
			  const int numBasis,
			  const int numQuadPts);

  /** Compute total strain in at quadrature points of a cell.
   *
   * @param strain Strain tensor at quadrature points.
   * @param basisDeriv Derivatives of basis functions at quadrature points.
   * @param disp Displacement at vertices of cell.
   * @param numBasis Number of basis functions for cell.
   * @param numQuadPts Number of quadrature points.
   */
  static
  void _calcTotalStrain3D(double_array* strain,
			  const double_array& basisDeriv,
			  const double_array& disp,
			  const int numBasis,
			  const int numQuadPts);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  /// Elastic material associated with integrator
  materials::ElasticMaterial* _material;

  /// Buffers for output.
  topology::Fields<topology::Field<topology::Mesh> >* _outputFields;

// PROTECTED TYPEDEFS ///////////////////////////////////////////////////
protected :

  typedef pylith::topology::Mesh::SieveMesh SieveMesh;
  typedef pylith::topology::Mesh::RealSection RealSection;

  typedef pylith::topology::Field<pylith::topology::Mesh>::RestrictVisitor RestrictVisitor;
  typedef pylith::topology::Field<pylith::topology::Mesh>::UpdateAddVisitor UpdateAddVisitor;
  typedef ALE::ISieveVisitor::IndicesVisitor<RealSection,SieveMesh::order_type,PetscInt> IndicesVisitor;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented.
  IntegratorElasticity(const IntegratorElasticity&);

  /// Not implemented
  const IntegratorElasticity& operator=(const IntegratorElasticity&);

}; // IntegratorElasticity

#endif // pylith_feassemble_integratorelasticity_hh


// End of file 
