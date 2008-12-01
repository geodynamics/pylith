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

#include "FaultCohesive.hh" // ISA FaultCohesive
#include "pylith/feassemble/Integrator.hh" // ISA Integrator

#include <map> // HASA std::map
#include <string> // HASA std::string

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class FaultCohesiveKin;
    class TestFaultCohesiveKin; // unit testing

    class EqKinSrc; // HOLDSA EqKinSrc
  } // faults
} // pylith

/*
namespace spatialdata {
  namespace geocoords {
    class CoordSys; // USES CoordSys
  } // geocoords
} // spatialdata
*/


/// C++ implementation for a fault surface with kinematic (prescribed)
/// slip implemented with cohesive elements.
class pylith::faults::FaultCohesiveKin : public FaultCohesive,
					 public feassemble::Integrator
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
   * @param sources Array of kinematic earthquake sources.
   * @param numSources Number of earthquake sources
   */
  void eqsrcs(const char** names,
	      EqKinSrc** sources,
	      const int numSources);

  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh PETSc mesh
   * @param cs Coordinate system for mesh
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction; only applies to fault surfaces in a 3-D domain).
   * @param normalDir General preferred direction for fault normal
   *   (used to pick which of two possible normal directions for
   *   interface; only applies to fault surfaces in a 3-D domain).
   * @param matDB Database of bulk elastic properties for fault region
   *   (used to improve conditioning of Jacobian matrix)
   */
  void initialize(const ALE::Obj<Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const double_array& upDir,
		  const double_array& normalDir,
		  spatialdata::spatialdb::SpatialDB* matDB);

  /** Integrate contributions to residual term (r) for operator that
   * require assembly across processors.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void integrateResidual(const ALE::Obj<real_section_type>& residual,
			 const double t,
			 topology::FieldsManager* const fields,
			 const ALE::Obj<Mesh>& mesh,
			 const spatialdata::geocoords::CoordSys* cs);

  /** Integrate contributions to residual term (r) for operator that
   * do not require assembly across cells, vertices, or processors.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void integrateResidualAssembled(const ALE::Obj<real_section_type>& residual,
				  const double t,
				  topology::FieldsManager* const fields,
				  const ALE::Obj<Mesh>& mesh,
				  const spatialdata::geocoords::CoordSys* cs);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator that do not require assembly across cells, vertices, or
   * processors.
   *
   * @param mat Sparse matrix
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void integrateJacobianAssembled(PetscMat* mat,
				  const double t,
				  topology::FieldsManager* const fields,
				  const ALE::Obj<Mesh>& mesh);

  /** Update state variables as needed.
   *
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void updateState(const double t,
		   topology::FieldsManager* const fields,
		   const ALE::Obj<Mesh>& mesh);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const ALE::Obj<Mesh>& mesh) const;

  /** Get vertex field associated with integrator.
   *
   * @param fieldType Type of field.
   * @param name Name of vertex field.
   * @param mesh PETSc mesh for problem. 
   * @param fields Fields manager.
   * @returns Vertex field.
   */
  const ALE::Obj<real_section_type>&
  vertexField(VectorFieldEnum* fieldType,
	      const char* name,
	      const ALE::Obj<Mesh>& mesh,
	      topology::FieldsManager* const fields);

  /** Get cell field associated with integrator.
   *
   * @param fieldType Type of field.
   * @param name Name of vertex field.
   * @param mesh PETSc mesh for problem.
   * @param fields Fields manager.
   * @returns Cell field.
   */
  const ALE::Obj<real_section_type>&
  cellField(VectorFieldEnum* fieldType,
	    const char* name,
	    const ALE::Obj<Mesh>& mesh,
	    topology::FieldsManager* const fields);

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
  void _calcOrientation(const double_array& upDir,
			const double_array& normalDir);

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
  void _calcTractionsChange(ALE::Obj<real_section_type>* tractions,
			    const ALE::Obj<Mesh>& mesh,
			    const ALE::Obj<real_section_type>& solution);

  /// Allocate scalar field for output of vertex information.
  void _allocateBufferVertexScalar(void);

  /// Allocate vector field for output of vertex information.
  void _allocateBufferVertexVector(void);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesiveKin(const FaultCohesiveKin& m);

  /// Not implemented
  const FaultCohesiveKin& operator=(const FaultCohesiveKin& m);

  // PRIVATE TYPEDEFS ///////////////////////////////////////////////////
private :

  typedef std::map<std::string, EqKinSrc*> srcs_type;

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  srcs_type _eqSrcs; ///< Array of kinematic earthquake sources.
  
  /// Field over fault mesh vertices of pseudo-stiffness values for
  /// scaling constraint information to improve conditioning of
  /// Jacobian matrix.
  ALE::Obj<real_section_type> _pseudoStiffness;

  /// Field over fault mesh vertices of area associated with each vertex.
  ALE::Obj<real_section_type> _area;

  /// Field over the fault mesh vertices of orientation of fault
  /// surface.
  ALE::Obj<real_section_type> _orientation;

  /// Field over the fault mesh vertices of vector field of current
  /// slip or slip increment.
  ALE::Obj<real_section_type> _slip;

  /// Field over the fault mesh vertices of vector field of cumulative slip.
  ALE::Obj<real_section_type> _cumSlip;

  /// Map label of cohesive cell to label of cells in fault mesh.
  std::map<Mesh::point_type, Mesh::point_type> _cohesiveToFault;

  /// Scalar field for vertex information over fault mesh.
  ALE::Obj<real_section_type> _bufferVertexScalar;

  /// Vector field for vertex information over fault mesh.
  ALE::Obj<real_section_type> _bufferVertexVector;

  /// Handle to field managed elsewhere for data over fault mesh.
  ALE::Obj<real_section_type> _bufferTmp;

}; // class FaultCohesiveKin

#include "FaultCohesiveKin.icc" // inline methods

#endif // pylith_faults_faultcohesivekin_hh


// End of file 
