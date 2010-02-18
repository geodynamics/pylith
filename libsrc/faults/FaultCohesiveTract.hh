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

/** @file libsrc/faults/FaultCohesiveTract.hh
 *
 * @brief C++ implementation for a fault surface with tractions
 * applied to the fault surface using cohesive cells.
 */

#if !defined(pylith_faults_faultcohesivetract_hh)
#define pylith_faults_faultcohesivetract_hh

// Include directives ---------------------------------------------------
#include "FaultCohesive.hh" // ISA FaultCohesive

// FaultCohesiveTract -----------------------------------------------------
/** 
 * @brief C++ implementation for a fault surface with tractions
 * applied to the fault surface using cohesive cells.
 */
class pylith::faults::FaultCohesiveTract : public FaultCohesive
{ // class FaultCohesiveTract
  friend class TestFaultCohesiveTract; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesiveTract(void);

  /// Destructor.
  virtual
  ~FaultCohesiveTract(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);

  /** Sets the spatial database for the inital tractions
   * @param dbs spatial database for initial tractions
   */
  void dbInitial(spatialdata::spatialdb::SpatialDB* dbs);
  
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
   */
  void initialize(const topology::Mesh& mesh,
		  const double upDir[3],
		  const double normalDir[3]);

  /** Integrate contribution of cohesive cells to residual term.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidual(const topology::Field<topology::Mesh>& residual,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateJacobian(topology::Jacobian* jacobian,
			 const double t,
			 topology::SolutionFields* const fields);
  
  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const topology::Mesh& mesh) const;

  /** Get vertex field associated with integrator.
   *
   * @param name Name of vertex field.
   * @param fields Solution fields.
   *
   * @returns Vertex field.
   */
  const topology::Field<topology::SubMesh>&
  vertexField(const char* name,
	      const topology::SolutionFields* fields =0);
  
  /** Get cell field associated with integrator.
   *
   * @param name Name of cell field.
   * @param fields Solution fields.
   *
   * @returns Cell field.
   */
  const topology::Field<topology::SubMesh>&
  cellField(const char* name,
	    const topology::SolutionFields* fields =0);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Calculate orientation at quadrature points.
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

  /** Get initial tractions using a spatial database.
   */
  void _getInitialTractions(void);

  /** Setup fault constitutive model.
   */
  void _initConstitutiveModel(void);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Database for initial tractions.
  spatialdata::spatialdb::SpatialDB* _dbInitial;


  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesiveTract(const FaultCohesiveTract&);

  /// Not implemented
  const FaultCohesiveTract& operator=(const FaultCohesiveTract&);

}; // class FaultCohesiveTract

#endif // pylith_faults_faultcohesivetract_hh


// End of file 
