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

/** @file libsrc/materials/EffectiveStress.hh
 *
 * @brief C++ EffectiveStress object.
 *
 * This class contains bracketing and root-finding functions for
 * materials that use an effective stress formulation.
 */

#if !defined(pylith_materials_effectivestress_hh)
#define pylith_materials_effectivestress_hh

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class EffectiveStress;
  } // materials

} // pylith

/// C++ abstract base class for Material object.
class pylith::materials::EffectiveStress
{ // class EffectiveStress

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Get effective stress from initial guess.
   *
   * The stressScale argument should provide a reasonable initial
   * guess in the case where the
   * actual initial guess is zero.
   *
   * @param effStressInitialGuess Initial guess for effective stress.
   * @param effStressParams Parameters used in computing effective stress.
   * @param material Material with effective stress function.
   *
   * @returns Computed effective stress.
   */
  template<typename material_type>
  static
  double getEffStress(const double effStressInitialGuess,
		      const double stressScale,
		      material_type* const material);

  // PRIVATE METHODS /////////////////////////////////////////////////////
private :

  /** Bracket effective stress.
   *
   * @param px1 Initial guess for first bracket.
   * @param px2 Initial guess for second bracket.
   * @param material Material with effective stress function.
   *
   */
  template<typename material_type>
  static
  void _bracketEffStress(double* px1,
			 double* px2,
			 material_type* const material);

  /** Solve for effective stress using Newton's method with bisection.
   *
   * @param x1 Initial guess for first bracket.
   * @param x2 Initial guess for second bracket.
   * @param material Material with effective stress function.
   *
   * @returns Computed effective stress.
   */
  template<typename material_type>
  static
  double _findEffStress(double x1,
			double x2,
			material_type* const material);

}; // class EffectiveStress

#endif // pylith_materials_effectivestress_hh

#include "EffectiveStress.icc" // template methods

// End of file 
