// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/FaultCohesiveLagrange.hh
 *
 * @brief C++ abstract base class for implementing fault slip using
 * Lagrange multipliers and cohesive cells.
 */

#if !defined(pylith_faults_faultcohesivelagrange_hh)
#define pylith_faults_faultcohesivelagrange_hh

// Include directives ---------------------------------------------------
#include "FaultCohesive.hh" // ISA FaultCohesive

// FaultCohesiveLagrange -----------------------------------------------------
/**
 * @brief C++ abstract base class for implementing falt slip using
 * Lagrange multipliers and cohesive cells.
 *
 * Fault boundary condition is specified using Lagrange
 * multipliers. The constraints are associated with "constraint"
 * vertices which sit between the pair of vertices on each side of the
 * fault. 
 *
 * The ordering of vertices in a cohesive cell is the vertices on the
 * "negative" side of the fault, the corresponding entries on the
 * "positive" side of the fault, and then the corresponding constraint
 * vertices.
 *
 * The system without Lagrange multipliers is
 *
 * [A(t+dt)]{u(t+dt)} = {b(t+dt)}
 *
 * With Lagrange multipliers this system becomes
 *
 * [A(t+dt) C^T ]{ u(t+dt) } = {b(t+dt)}
 * [ C      0   ]{ L(t+dt) }   {D(t+dt)}
 *
 * where C is the matrix of Lagrange constraints, L is the vector of
 * Lagrange multiplies (internal forces in this case), and D is the
 * fault slip.
 *
 * We solve for the increment in the displacement field, so we rewrite
 * the system as
 *
 * [A(t+dt) C^T ]{ du(t) } = {b(t+dt)} - [A(t+dt) C^T ]{ u(t) }
 * [ C      0   ]{ dL(t) }   {D(t+dt)}   [ C      0   ]{ L(t) }
 * 
 * We form the residual as
 *
 * {r(t+dt)} = {b(t+dt)} - [A(t+dt) C^T ]{ u(t)+du(t) }
 *             {D(t+dt)}   [ C      0   ]{ L(t)+dL(t) }
 * 
 * The terms in the residual contributing to the DOF at the Lagrange
 * vertices are
 *
 * {r(t+dt)} = {D(t+dt)} - [C]{u(t)+dt(t)}
 *
 * The term in the residual contributing to the DOF at the
 * non-Lagrange vertices of the cohesive cells is
 *
 * {r(t+dt)} = -[C]^T{L(t)+dL(t)}
 *
 * We integrate the Lagrange multiplier term and the relative
 * displacement term over the cohesive cells, because this introduces
 * weighting of the orientation of the fault for the direction of slip
 * at the vertices of the cohesive cells.
 */
class pylith::faults::FaultCohesiveLagrange : public FaultCohesive
{ // class FaultCohesiveLagrange
  friend class TestFaultCohesiveLagrange; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesiveLagrange(void);

  /// Destructor.
  virtual
  ~FaultCohesiveLagrange(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh Finite-element mesh.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction; applies to fault surfaces in 2-D and 3-D).
   */
  virtual
  void initialize(const topology::Mesh& mesh,
		  const PylithScalar upDir[3]);

  /** Split solution field for separate preconditioning.
   *
   * @param field Solution field.
   */
  void splitField(topology::Field<topology::Mesh>* field);

  /** Integrate contributions to residual term (r) for operator that
   * require assembly processors.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  virtual
  void integrateResidual(const topology::Field<topology::Mesh>& residual,
				  const PylithScalar t,
				  topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator that require assembly across processors.
   *
   * @param jacobian Sparse matrix
   * @param t Current time
   * @param fields Solution fields
   */
  virtual
  void integrateJacobian(topology::Jacobian* jacobian,
				  const PylithScalar t,
				  topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator that require assembly processors.
   *
   * @param jacobian Diagonal Jacobian matrix as a field.
   * @param t Current time
   * @param fields Solution fields
   */
  virtual
  void integrateJacobian(topology::Field<topology::Mesh>* jacobian,
				  const PylithScalar t,
				  topology::SolutionFields* const fields);

  /** Compute custom fault precoditioner using Schur complement.
   *
   * We have J = [A C^T]
   *             [C   0]
   *
   * We approximate C A^(-1) C^T.
   *
   * @param pc PETSc preconditioner structure.
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param fields Solution fields
   */
  virtual
  void calcPreconditioner(PetscMat* const precondMatrix,
			  topology::Jacobian* const jacobian,
			  topology::SolutionFields* const fields);

  /** Adjust solution from solver with lumped Jacobian to match Lagrange
   *  multiplier constraints.
   *
   * @param fields Solution fields
   * @param jacobian Jacobian of the system.
   */
  virtual
  void adjustSolnLumped(topology::SolutionFields* fields,
			const topology::Field<topology::Mesh>& jacobian);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  virtual
  void verifyConfiguration(const topology::Mesh& mesh) const;

  /** Verify constraints are acceptable.
   *
   * @param field Solution field.
   */
  virtual
  void checkConstraints(const topology::Field<topology::Mesh>& solution) const;

  /** Get cell field associated with integrator.
   *
   * @param name Name of cell field.
   * @param fields Solution fields.
   * @returns Cell field.
   */
  const topology::Field<topology::SubMesh>&
  cellField(const char* name,
	    const topology::SolutionFields* fields =0);

  // PROTECTED STRUCTS //////////////////////////////////////////////////
protected :

  /** Data structure to hold relations between vertices in cohesive cell
   *  in mesh of domain and cell of fault mesh.
   */
  struct CohesiveInfo {
    int lagrange; ///< Vertex associated with Lagrange multiplier.
    int positive; ///< Vertex on positive side of the fault.
    int negative; ///< Vertex on negative side of the fault.
    int fault; ///< Vertex in fault mesh.
  };

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Initialize auxiliary cohesive cell information.
   *
   * @param mesh Finite-element mesh of the domain.
   */
  void _initializeCohesiveInfo(const topology::Mesh& mesh);

  /** Compute change in tractions on fault surface using solution.
   *
   * @param tractions Field for tractions.
   * @param solution Solution over domain
   */
  void _calcTractionsChange(topology::Field<topology::SubMesh>* tractions,
          const topology::Field<topology::Mesh>& solution);

  /** Transform field from local (fault) coordinate system to
   * global coordinate system.
   *
   * @param field Field to transform.
   */
  void _faultToGlobal(topology::Field<topology::SubMesh>* field);

  /** Transform field from global coordinate system to local (fault)
   * coordinate system.
   *
   * @param field Field to transform.
   */
  void _globalToFault(topology::Field<topology::SubMesh>* field);

  /// Allocate buffer for vector field.
  void _allocateBufferVectorField(void);

  /// Allocate buffer for scalar field.
  void _allocateBufferScalarField(void);

  /** Get submatrix of Jacobian matrix associated with the negative
   *  and positive sides of the fault.
   *
   *  The submatrix and the mapping from global indices to submatrix
   *  indices are returned through the arguments.
   *
   * @param submatrix PETSc matrix for submatrix.
   * @param mapGlobalToLocal Mapping from global to local indices.
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param fields Solution fields
   */
  void _getJacobianSubmatrixNP(PetscMat* submatrix,
			       std::map<int,int>* mapGlobalToLocal,
			       const topology::Jacobian& jacobian,
			       const topology::SolutionFields& fields);

  /** Adjust solution in lumped formulation to match slip for 1-D.
   *
   * @param lagrangeIncr Increment in Lagrange multiplier.
   * @param dispTIncrN Disp increment assoc. w/vertex on - side of fault.
   * @param dispTIncrP Disp increment assoc. w/vertex on + side of fault.
   * @param slip Slip associated with Lagrange multipler vertex.
   * @param orientation Orientation associated with Lagrange multipluer vertex.
   * @param dispTN Displacement associated with vertex on - side of fault.
   * @param dispTP Displacement associated with vertex on + side of fault.
   * @param residualN Residual associated with vertex on - side of fault.
   * @param residualP Residual associated with vertex on + side of fault.
   * @param jacobianN Jacobian associated with vertex on - side of fault.
   * @param jacobianP Jacobian associated with vertex on + side of fault.
   */
  void _adjustSolnLumped1D(scalar_array* lagrangeIncr,
			   scalar_array* dispTIncrN,
			   scalar_array* dispTIncrP,
			   const scalar_array& slip,
			   const scalar_array& orientation,
			   const scalar_array& dispTN,
			   const scalar_array& dispTP,
			   const scalar_array& residualN,
			   const scalar_array& residualP,
			   const scalar_array& jacobianN,
			   const scalar_array& jacobianP);

  /** Adjust solution in lumped formulation to match slip for 2-D.
   *
   * @param lagrangeIncr Increment in Lagrange multiplier.
   * @param dispTIncrN Disp increment assoc. w/vertex on - side of fault.
   * @param dispTIncrP Disp increment assoc. w/vertex on + side of fault.
   * @param slip Slip associated with Lagrange multipler vertex.
   * @param orientation Orientation associated with Lagrange multipluer vertex.
   * @param dispTN Displacement associated with vertex on - side of fault.
   * @param dispTP Displacement associated with vertex on + side of fault.
   * @param residualN Residual associated with vertex on - side of fault.
   * @param residualP Residual associated with vertex on + side of fault.
   * @param jacobianN Jacobian associated with vertex on - side of fault.
   * @param jacobianP Jacobian associated with vertex on + side of fault.
   */
  void _adjustSolnLumped2D(scalar_array* lagrangeIncr,
			   scalar_array* dispTIncrN,
			   scalar_array* dispTIncrP,
			   const scalar_array& slip,
			   const scalar_array& orientation,
			   const scalar_array& dispTN,
			   const scalar_array& dispTP,
			   const scalar_array& residualN,
			   const scalar_array& residualP,
			   const scalar_array& jacobianN,
			   const scalar_array& jacobianP);

  /** Adjust solution in lumped formulation to match slip for 3-D.
   *
   * @param lagrangeIncr Increment in Lagrange multiplier.
   * @param dispTIncrN Disp increment assoc. w/vertex on - side of fault.
   * @param dispTIncrP Disp increment assoc. w/vertex on + side of fault.
   * @param slip Slip associated with Lagrange multipler vertex.
   * @param orientation Orientation associated with Lagrange multipluer vertex.
   * @param dispTN Displacement associated with vertex on - side of fault.
   * @param dispTP Displacement associated with vertex on + side of fault.
   * @param residualN Residual associated with vertex on - side of fault.
   * @param residualP Residual associated with vertex on + side of fault.
   * @param jacobianN Jacobian associated with vertex on - side of fault.
   * @param jacobianP Jacobian associated with vertex on + side of fault.
   */
  void _adjustSolnLumped3D(scalar_array* lagrangeIncr,
			   scalar_array* dispTIncrN,
			   scalar_array* dispTIncrP,
			   const scalar_array& slip,
			   const scalar_array& orientation,
			   const scalar_array& dispTN,
			   const scalar_array& dispTP,
			   const scalar_array& residualN,
			   const scalar_array& residualP,
			   const scalar_array& jacobianN,
			   const scalar_array& jacobianP);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  /// Array of cohesive vertex information.
  std::vector<CohesiveInfo> _cohesiveVertices;

  /// Map label of cohesive cell to label of cells in fault mesh.
  std::map<topology::Mesh::SieveMesh::point_type,
           topology::SubMesh::SieveMesh::point_type> _cohesiveToFault;

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /// Initialize logger.
  void _initializeLogger(void);
  
  /** Calculate orientation at fault vertices.
   *
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction; applies to fault surfaces in 2-D and 3-D).
   */
  void _calcOrientation(const PylithScalar upDir[3]);

  /// Calculate fault area field.
  void _calcArea(void);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesiveLagrange(const FaultCohesiveLagrange&);

  /// Not implemented
  const FaultCohesiveLagrange& operator=(const FaultCohesiveLagrange&);

}; // class FaultCohesiveLagrange

#endif // pylith_faults_faultcohesivelagrange_hh


// End of file 
