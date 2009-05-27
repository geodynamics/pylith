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

  // PUBLIC TYPEDEFS ///////////////////////////////////////////////////
public :

  /// Member prototype for effStressFunc()
  typedef double (*effStressFunc_fn_type)
  (const double,
   const EffStressStruct&);

  /// Member prototype for effStressDFunc()
  typedef double (*effStressDFunc_fn_type)
  (const double,
   const EffStressStruct&);
  
  /// Member prototype for effStressFuncDFunc()
  typedef void (*effStressFuncDFunc_fn_type)
  (const double,
   const EffStressStruct&,
   double*,
   double*);

  // PUBLIC MEMBERS ///////////////////////////////////////////////////
public :

  /// Metod to use for effStressFunc().
  effStressFunc_fn_type effStressFunc;

  /// Metod to use for effStressDFunc().
  effStressDFunc_fn_type effStressDFunc;

  /// Metod to use for effStressFuncDFunc().
  effStressFuncDFunc_fn_type effStressFuncDFunc;


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
  static
  double getEffStress(const double effStressInitialGuess,
		      const double stressScale,
		      const EffStressStruct& effStressParams,
		      effStressFunc_fn_type effStressFunc,
		      effStressFuncDFunc_fn_type effStressFuncDFunc);

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
  static
  void _bracketEffStress(double* px1,
			 double* px2,
			 const EffStressStruct& effStressParams,
			 effStressFunc_fn_type effStressFunc);

  /** Solve for effective stress using Newton's method with bisection.
   *
   * @param x1 Initial guess for first bracket.
   * @param x2 Initial guess for second bracket.
   * @param effStressParams Parameters used in computing effective stress.
   * @param effStressFunc Function to compute effective stress only.
   * @param effStressFuncDFunc Function to compute effective stress and derivative.
   *
   * @returns Computed effective stress.
   */
  static
  double _findEffStress(double x1,
			double x2,
			const EffStressStruct& effStressParams,
			effStressFunc_fn_type effStressFunc,
			effStressFuncDFunc_fn_type effStressFuncDFunc);

}; // class EffectiveStress

#endif // pylith_materials_effectivestress_hh


// End of file 
