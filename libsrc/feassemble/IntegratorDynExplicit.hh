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
 * @file pylith/feassemble/IntegratorDynExplicit.hh
 *
 * @brief Abstract base class for explicit time integration of
 * finite-element actions.
 *
 * Computes terms A and b in A(t) u(t+dt) = b(u(t), u(t-dt)), where
 * A(t) is a sparse matrix or vector, u(t+dt) is the field we want to
 * compute at time t+dt and b is a vector that depends on the field at
 * time t and t-dt.
 */

#if !defined(pylith_feassemble_integratordynexplicit_hh)
#define pylith_feassemble_integratordynexplicit_hh

#include <petscmesh.h> // USES Mesh
#include "pylith/utils/petscfwd.h" // USES PetscMat

namespace pylith {
  namespace feassemble {
    class IntegratorDynExplicit;
    class TestIntegratorDynExplicit;

    class Quadrature; // HOLDSA Quadrature
  } // feassemble
} // pylith

namespace spatialdata {
  namespace spatialdb {
    class SpatialDB; // USES SpatialDB
  } // spatialdb
  namespace geocoords {
    class CoordSys; // USES CoordSys
  } // geocoords
} // spatialdata

class pylith::feassemble::IntegratorDynExplicit
{ // Integrator
  friend class TestIntegratorDynExplicit; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :

  typedef ALE::Mesh Mesh;
  typedef Mesh::topology_type topology_type;
  typedef topology_type::point_type point_type;
  typedef Mesh::real_section_type real_section_type;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  IntegratorDynExplicit(void);

  /// Destructor
  virtual
  ~IntegratorDynExplicit(void);

  /// Create a copy of this object.
  virtual
  IntegratorDynExplicit* clone(void) const = 0;

  /** Integrate residual term (b) for dynamic elasticity term 
   * for 3-D finite elements.
   *
   * @param fieldOut Output field
   * @param fieldInT Input field at time t
   * @param fieldInTmdt Input field at time t-dt
   * @param coordinates Field of cell vertex coordinates
   */
  virtual 
  void integrateResidual(const ALE::Obj<real_section_type>& fieldOut,
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
   * @param fieldIn Input field at time t
   * @param coordinates Field of cell vertex coordinates
   */
  virtual 
  void integrateJacobian(const ALE::Obj<real_section_type>& fieldOut,
			 const ALE::Obj<real_section_type>& fieldIn,
			 const ALE::Obj<real_section_type>& coordinates) = 0;
  
  /** Set quadrature for integrating finite-element quantities.
   *
   * @param q Quadrature for integrating.
   */
  void quadrature(const Quadrature* q);

  /** Initialize, get material property parameters from database.
   *
   * @param mesh PETSc mesh
   * @param cs Pointer to coordinate system of vertices
   * @param db Pointer to spatial database with material property parameters
   */
  virtual
  void initialize(ALE::Obj<ALE::Mesh>& mesh,
		  spatialdata::geocoords::CoordSys* cs,
		  spatialdata::spatialdb::SpatialDB* db) = 0;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param i Integrator to copy
   */
  IntegratorDynExplicit(const IntegratorDynExplicit& i);

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

  /// Not implemented
  const IntegratorDynExplicit& operator=(const IntegratorDynExplicit&);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  Quadrature* _quadrature; ///< Quadrature for integrating finite-element

  /// Vector local to cell containing result of integration action
  real_section_type::value_type* _cellVector;

  /// Matrix local to cell containing result of integration
  real_section_type::value_type* _cellMatrix;

}; // IntegratorDynExplicit

#endif // pylith_feassemble_integratordynexplicit_hh


// End of file 
