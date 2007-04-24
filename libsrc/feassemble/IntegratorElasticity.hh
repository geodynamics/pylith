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
 * @file pylith/feassemble/IntegratorElasticity.hh
 *
 * @brief C++ Utility class for general operations associated with
 * integrating finite-element actions for solid elasticity.
 */

#if !defined(pylith_feassemble_integratorelasticity_hh)
#define pylith_feassemble_integratorelasticity_hh

#include "pylith/utils/array.hh" // USES double_array, std::vector

namespace pylith {
  namespace feassemble {
    class IntegratorElasticity;
    class TestIntegratorElasticity;
  } // feassemble
} // pylith

class pylith::feassemble::IntegratorElasticity
{ // IntegratorElasticity
  friend class TestIntegratorElasticity; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Compute total strain in at quadrature points of a cell.
   *
   * @param strain Strain tensor at quadrature points.
   * @param basisDeriv Derivatives of basis functions at quadrature points.
   * @param disp Displacement at vertices of cell.
   * @param dimension Dimension of cell.
   * @param numBasis Number of basis functions for cell.
   */

  static
  void calcTotalStrain1D(std::vector<double_array>* strain,
			 const double_array& basisDeriv,
			 const double* disp,
			 const int numBasis);

  /** Compute total strain in at quadrature points of a cell.
   *
   * @param strain Strain tensor at quadrature points.
   * @param basisDeriv Derivatives of basis functions at quadrature points.
   * @param disp Displacement at vertices of cell.
   * @param numBasis Number of basis functions for cell.
   */

  static
  void calcTotalStrain2D(std::vector<double_array>* strain,
			 const double_array& basisDeriv,
			 const double* disp,
			 const int numBasis);

  /** Compute total strain in at quadrature points of a cell.
   *
   * @param strain Strain tensor at quadrature points.
   * @param basisDeriv Derivatives of basis functions at quadrature points.
   * @param disp Displacement at vertices of cell.
   * @param numBasis Number of basis functions for cell.
   */

  static
  void calcTotalStrain3D(std::vector<double_array>* strain,
			 const double_array& basisDeriv,
			 const double* disp,
			 const int numBasis);

}; // IntegratorElasticity

#endif // pylith_feassemble_integratorelasticity_hh

// End of file 
