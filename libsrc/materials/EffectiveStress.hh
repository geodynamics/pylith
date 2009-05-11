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
 * This class contains bracketing and root-finding functions for materials that
 * use an effective stress formulation.
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

  // PUBLIC STRUCTS /////////////////////////////////////////////////////
public :

  struct EffStressStruct {
    double stressScale;
    double ae;
    double b;
    double c;
    double d;
    double alpha;
    double dt;
    double effStressT;
    double powerLawExp;
    double viscosityCoeff;
  };

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Get effective stress from initial guess.
   *
   * @param effStressInitialGuess Initial guess for effective stress.
   * @param effStressParams Parameters used in computing effective stress.
   * @param effStressFunc Function to compute effective stress only.
   * @param effStressFuncDFunc Function to compute effective stress and derivative.
   *
   * @returns Computed effective stress.
   */
  static double getEffStress(const double effStressInitialGuess,
			     const EffStressStruct& effStressParams,
			     effStressFuncType* effStressFunc,
			     effStressFuncDFuncType* effStressFuncDFunc);

  // PRIVATE METHODS /////////////////////////////////////////////////////
private :

  /** Bracket effective stress.
   *
   * @param x1 Initial guess for first bracket.
   * @param x2 Initial guess for second bracket.
   * @param effStressParams Parameters used in computing effective stress.
   * @param effStressFunc Function to compute effective stress only.
   *
   */
  void bracketEffStress(double x1,
			double x2,
			const double effStressParams,
			static double &effStressFunc);

  /** Solve for effective stress using Newton's method with bisection.
   *
   * @param x1 Initial guess for first bracket.
   * @param x2 Initial guess for second bracket.
   * @param effStressParams Parameters used in computing effective stress.
   * @param effStressFunc Function to compute effective stress only.
   *
   */
  void bracketEffStress(double x1,
			double x2,
			const double effStressParams,
			static double &effStressFunc);
}; // class EffectiveStress

#endif // pylith_materials_effectivestress_hh


// End of file 
