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

/** @file libsrc/bc/BoundaryCondition.hh
 *
 * @brief C++ abstract base class for BoundaryCondition object.
 *
 * Interface definition for boundary conditions.
 */

#if !defined(pylith_bc_boundarycondition_hh)
#define pylith_bc_boundarycondition_hh

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh, real_section_type
#include "pylith/utils/petscfwd.h" // USES PETScMat

#include <string> // HASA std::string

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class BoundaryCondition;
    class TestBoundaryCondition; // unit testing
  } // bc
} // pylith

/// Namespace for spatialdata package
namespace spatialdata {
  namespace geocoords {
    class CoordSys;
  } // geocoords

  namespace spatialdb {
    class SpatialDB;
  } // spatialdb
} // spatialdata

/// C++ abstract base class for BoundaryCondition object.
class pylith::bc::BoundaryCondition
{ // class BoundaryCondition
  friend class TestBoundaryCondition; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  BoundaryCondition(void);

  /// Destructor.
  virtual
  ~BoundaryCondition(void);

  /** Set identifier of fault.
   *
   * @param value BoundaryCondition identifier
   */
  void id(const int value);

  /** Get identifier of fault.
   *
   * @returns BoundaryCondition identifier
   */
  int id(void) const;

  /** Set label of fault.
   *
   * @param value Label of fault
   */
  void label(const char* value);

  /** Get label of fault.
   *
   * @returns Label of fault
   */
  const std::string& label(void) const;

  /** Set database for boundary condition parameters.
   *
   * @param db Spatial database
   */
  void db(spatialdata::spatialdb::SpatialDB* const db);

  /** Initialize boundary condition.
   *
   * @param mesh PETSc mesh
   * @param cs Coordinate system for mesh
   */
  virtual
  void initialize(const ALE::Obj<ALE::Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs) = 0;

  /** Set number of degrees of freedom that are constrained at points in field.
   *
   * @param field Solution field
   * @param mesh PETSc mesh
   */
  virtual
  void setConstraintSizes(const ALE::Obj<real_section_type>& field,
			  const ALE::Obj<ALE::Mesh>& mesh);

  /** Set which degrees of freedom are constrained at points in field.
   *
   * @param field Solution field
   * @param mesh PETSc mesh
   */
  virtual
  void setConstraints(const ALE::Obj<real_section_type>& field,
		      const ALE::Obj<ALE::Mesh>& mesh);

  /** Set constrained degrees of freedom in field.
   *
   * @param jacobian Matrix for Jacobian of system
   * @param field Solution field
   * @param mesh PETSc mesh
   */
  virtual
  void integrateJacobian(PetscMat* jacobian,
			 const ALE::Obj<real_section_type>& field,
			 const ALE::Obj<ALE::Mesh>& mesh);

  /** Set constrained degrees of freedom in field.
   *
   * @param residual Residual field
   * @param fieldT Solution field at time t
   * @param fieldTmdt Solution field at time t-dt
   * @param mesh PETSc mesh
   */
  virtual
  void integrateResidual(const ALE::Obj<real_section_type>& residual,
			 const ALE::Obj<real_section_type>& fieldT,
			 const ALE::Obj<real_section_type>& fieldTmdt,
			 const ALE::Obj<ALE::Mesh>& mesh);

  /** Set values in field.
   *
   * @param t Current time
   * @param field Solution field
   * @param mesh PETSc mesh
   */
  virtual
  void setField(const double t,
		const ALE::Obj<real_section_type>& field,
		const ALE::Obj<ALE::Mesh>& mesh);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  BoundaryCondition(const BoundaryCondition& m);

  /// Not implemented
  const BoundaryCondition& operator=(const BoundaryCondition& m);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  int _id; ///< BoundaryCondition identifier
  std::string _label; ///< Label of fault
  spatialdata::spatialdb::SpatialDB* _db; ///< Spatial database w/parameters

}; // class BoundaryCondition

#include "BoundaryCondition.icc" // inline methods

#endif // pylith_bc_boundarycondition_hh


// End of file 
