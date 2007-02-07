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
 * Note: This object operates on a single finite-element family, which
 * is defined by the quadrature and a database of material property
 * parameters.
 *
 * Computes terms A and b in A(t) u(t+dt) = b(u(t), u(t-dt)), where
 * A(t) is a sparse matrix or vector, u(t+dt) is the field we want to
 * compute at time t+dt and b is a vector that depends on the field at
 * time t and t-dt.
 */

#if !defined(pylith_feassemble_integratordynexplicit_hh)
#define pylith_feassemble_integratordynexplicit_hh

#include "pylith/utils/petscfwd.h" // USES PetscMat

#include "Integrator.hh" // ISA Integrator

namespace pylith {
  namespace feassemble {
    class IntegratorDynExplicit;
    class TestIntegratorDynExplicit;
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

class pylith::feassemble::IntegratorDynExplicit : public Integrator
{ // Integrator
  friend class TestIntegratorDynExplicit; // unit testing

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
   * @param residual Residual field (output)
   * @param fieldInT Input field at time t
   * @param fieldInTmdt Input field at time t-dt
   * @param coordinates Field of cell vertex coordinates
   */
  virtual 
  void integrateResidual(const ALE::Obj<real_section_type>& residual,
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

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const IntegratorDynExplicit& operator=(const IntegratorDynExplicit&);

}; // IntegratorDynExplicit

#endif // pylith_feassemble_integratordynexplicit_hh


// End of file 
