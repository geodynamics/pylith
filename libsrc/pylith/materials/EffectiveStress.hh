// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/EffectiveStress.hh
 *
 * @brief C++ EffectiveStress object.
 */

#if !defined(pylith_materials_effectivestress_hh)
#define pylith_materials_effectivestress_hh

// Include directives ---------------------------------------------------
#include "materialsfwd.hh"

// EffectiveStress ------------------------------------------------------
/** @brief C++ EffectiveStress object.
 *
 * This class contains bracketing and root-finding functions for
 * materials that use an effective stress formulation.
 */
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
  PylithScalar calculate(const PylithScalar effStressInitialGuess,
		   const PylithScalar stressScale,
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
  void _bracket(PylithScalar* px1,
		PylithScalar* px2,
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
  PylithScalar _search(PylithScalar x1,
		 PylithScalar x2,
		 material_type* const material);

}; // class EffectiveStress

#endif // pylith_materials_effectivestress_hh

#include "EffectiveStress.icc" // template methods

// End of file 
