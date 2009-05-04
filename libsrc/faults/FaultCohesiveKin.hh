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
 * @brief C++ implementation for a fault surface with kinematic
 * (prescribed) slip implemented with cohesive elements.
 *
 * Fault boundary condition is specified using Lagrange
 * multipliers. The constraints are associated with "constraint"
 * vertices which sit between the pair of vertices on each side of the
 * fault. 
 *
 * The ordering of vertices in a cohesive cell is the vertices on the
 * one side of the fault, the corresponding entries on the other side
 * of the fault, and then the corresponding constraint vertices.
 *
 * [ K  aC^T ] [ U    ] = [ Fe ]
 * [ C   0   ] [ Fi/a ] = [ D  ]
 *
 * where K is the stiffness matrix, C is the matrix of Lagrange
 * constraints, U is the displacement field, Fe is the vector of
 * external forces, Fi is the vector of Lagrange multipers (forces), D
 * is the fault slip, and "a" is the conditioning value (taken to be
 * the shear modulus).
 */

#if !defined(pylith_faults_faultcohesivekin_hh)
#define pylith_faults_faultcohesivekin_hh

// Include directives ---------------------------------------------------
#include "FaultCohesive.hh" // ISA FaultCohesive

#include "pylith/topology/SubMesh.hh" // ISA Integrator<Quadrature<SubMesh> >
#include "pylith/feassemble/Quadrature.hh" // ISA Integrator<Quadrature>
#include "pylith/feassemble/Integrator.hh" // ISA Integrator

#include <map> // HASA std::map
#include <string> // HASA std::string

// FaultCohesiveKin -----------------------------------------------------
class pylith::faults::FaultCohesiveKin : public FaultCohesive,
					 public feassemble::Integrator<feassemble::Quadrature<topology::SubMesh> >
{ // class FaultCohesiveKin
  friend class TestFaultCohesiveKin; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesiveKin(void);

  /// Destructor.
  virtual
  ~FaultCohesiveKin(void);

  /** Set kinematic earthquake sources.
   *
   * @param names Array of kinematic earthquake source names.
   * @param numNames Number of earthquake sources.
   * @param sources Array of kinematic earthquake sources.
   * @param numSources Number of earthquake sources.
   */
  void eqsrcs(const char* const* names,
	      const int numNames,
	      EqKinSrc** sources,
	      const int numSources);

  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh Finite-element mesh.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction; only applies to fault surfaces in a 3-D domain).
   * @param normalDir General preferred direction for fault normal
   *   (used to pick which of two possible normal directions for
   *   interface; only applies to fault surfaces in a 3-D domain).
   * @param matDB Database of bulk elastic properties for fault region
   *   (used to improve conditioning of Jacobian matrix)
   */
  void initialize(const topology::Mesh& mesh,
		  const double upDir[3],
		  const double normalDir[3],
		  spatialdata::spatialdb::SpatialDB* matDB);

  /** Integrate contributions to residual term (r) for operator that
   * require assembly across processors.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidual(const topology::Field<topology::Mesh>& residual,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to residual term (r) for operator that
   * do not require assembly across cells, vertices, or processors.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidualAssembled(const topology::Field<topology::Mesh>& residual,
				  const double t,
				  topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator that do not require assembly across cells, vertices, or
   * processors.
   *
   * @param mat Sparse matrix
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void integrateJacobianAssembled(topology::Jacobian* jacobian,
				  const double t,
				  topology::SolutionFields* const fields);

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

  /** Get vertex field associated with integrator.
   *
   * @param name Name of cell field.
   * @param fields Solution fields.
   * @returns Vertex field.
   */
  const topology::Field<topology::SubMesh>&
  vertexField(const char* name,
	      const topology::SolutionFields* fields =0);

  /** Get cell field associated with integrator.
   *
   * @param name Name of cell field.
   * @param fields Solution fields.
   * @returns Cell field.
   */
  const topology::Field<topology::SubMesh>&
  cellField(const char* name,
	    const topology::SolutionFields* fields =0);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Cohesive cells use Lagrange multiplier constraints?
   *
   * @returns True if implementation using Lagrange multiplier
   * constraints, false otherwise.
   */
  bool _useLagrangeConstraints(void) const;

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Calculate orientation at fault vertices.
   *
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction; only applies to fault surfaces in a 3-D domain).
   * @param normalDir General preferred direction for fault normal
   *   (used to pick which of two possible normal directions for
   *   interface; only applies to fault surfaces in a 3-D domain).
   */
  void _calcOrientation(const double upDir[3],
			const double normalDir[3]);

  /** Calculate conditioning field.
   *
   * @param cs Coordinate system for mesh
   * @param matDB Database of bulk elastic properties for fault region
   *   (used to improve conditioning of Jacobian matrix)
   */
  void _calcConditioning(const spatialdata::geocoords::CoordSys* cs,
			 spatialdata::spatialdb::SpatialDB* matDB);

  /// Calculate fault area field.
  void _calcArea(void);

  /** Compute change in tractions on fault surface using solution.
   *
   * @param tractions Field for tractions.
   * @param mesh Finite-element mesh for domain
   * @param solution Solution over domain
   */
  void _calcTractionsChange(topology::Field<topology::SubMesh>* tractions,
			    const topology::Field<topology::Mesh>& solution);

  /// Allocate buffer for vector field.
  void _allocateBufferVectorField(void);

  /// Allocate buffer for scalar field.
  void _allocateBufferScalarField(void);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesiveKin(const FaultCohesiveKin&);

  /// Not implemented
  const FaultCohesiveKin& operator=(const FaultCohesiveKin&);

  // PRIVATE TYPEDEFS ///////////////////////////////////////////////////
private :

  typedef std::map<std::string, EqKinSrc*> srcs_type;

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  srcs_type _eqSrcs; ///< Array of kinematic earthquake sources.

  /// Fields for fault information.
  topology::Fields<topology::Field<topology::SubMesh> >* _fields;

  /// Buffer for vector field over fault vertices.
  topology::Field<topology::SubMesh>* _bufferVectorField;
  
  /// Buffer for scalar field over fault vertices.
  topology::Field<topology::SubMesh>* _bufferScalarField;
  
  /// Map label of cohesive cell to label of cells in fault mesh.
  std::map<topology::Mesh::SieveMesh::point_type, 
	   topology::SubMesh::SieveMesh::point_type> _cohesiveToFault;

}; // class FaultCohesiveKin

#include "FaultCohesiveKin.icc" // inline methods

#endif // pylith_faults_faultcohesivekin_hh


// End of file 
