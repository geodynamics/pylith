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
 * @brief Abstract base class for integration of finite-element
 * actions.
 *
 * Note: Each object operates on a single finite-element family, which
 * is defined by the quadrature and a database of material property
 * parameters.
 */

#if !defined(pylith_feassemble_integrator_hh)
#define pylith_feassemble_integrator_hh

#include <petscmesh.h> // USES Mesh

namespace pylith {
  namespace feassemble {
    class Integrator;
    class TestIntegrator;

    class Quadrature; // HOLDSA Quadrature
    class ParameterManager; // HOLDSA ParameterManager
  } // feassemble
} // pylith

class pylith::feassemble::Integrator
{ // Integrator
  friend class TestIntegrator; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :

  typedef ALE::Field::Mesh        Mesh;
  typedef Mesh::point_type        point_type;
  typedef Mesh::real_section_type real_section_type;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Integrator(void);

  /// Destructor
  virtual
  ~Integrator(void);

  /** Set quadrature for integrating finite-element
   * quantities. Quadrature should already be initialized.
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
