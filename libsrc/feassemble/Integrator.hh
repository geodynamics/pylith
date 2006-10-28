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
 * @file pylith/feassemble/Integrator.hh
 *
 * @brief Abstract base class for integration of finite-element actions.
 */

#if !defined(pylith_feassemble_integrator_hh)
#define pylith_feassemble_integrator_hh

#include <petscmesh.h> // USES Mesh
#include "pylith/utils/petscfwd.h" // USES PetscMat

namespace pylith {
  namespace feassemble {
    class Integrator;
    class TestIntegrator;

    class Quadrature; // HOLDSA Quadrature
  } // feassemble
} // pylith

class pylith::feassemble::Integrator
{ // Integrator
  friend class TestIntegrator; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :

  typedef ALE::Mesh Mesh;
  typedef Mesh::topology_type topology_type;
  typedef topology_type::point_type point_type;
  typedef Mesh::real_section_type real_section_type;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Integrator(void);

  /// Destructor
  virtual
  ~Integrator(void);

  /// Create a copy of this object.
  virtual
  Integrator* clone(void) const = 0;

  /** Integrate elasticity term for 3-D finite elements.
   *
   * @param fieldOut Output field
   * @param fieldIn Input field
   * @param coordinates Field of cell vertex coordinates
   */
  virtual 
  void integrateAction(const ALE::Obj<real_section_type>& fieldOut,
		       const ALE::Obj<real_section_type>& fieldIn,
		       const ALE::Obj<real_section_type>& coordinates) = 0;

  /** Compute matrix associated with operator.
   *
   * @param mat Sparse matrix
   * @param coordinates Field of cell vertex coordinates
   */
  virtual 
  void integrate(PetscMat* mat,
		 ALE::Obj<ALE::Mesh>& mesh,
		 const ALE::Obj<real_section_type>& fieldIn,
		 const ALE::Obj<real_section_type>& coordinates) = 0;

  /** Set quadrature for integrating finite-element quantities.
   *
   * @param q Quadrature for integrating.
   */
  void quadrature(const Quadrature* q);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param i Integrator to copy
   */
  Integrator(const Integrator& i);

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
  const Integrator& operator=(const Integrator&);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  Quadrature* _quadrature; ///< Quadrature for integrating finite-element

  /// Vector local to cell containing result of integration action
  real_section_type::value_type* _cellVector;

  /// Matrix local to cell containing result of integration
  real_section_type::value_type* _cellMatrix;

}; // Integrator

#endif // pylith_feassemble_integrator_hh

// End of file 
