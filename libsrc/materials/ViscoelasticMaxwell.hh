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

/** @file libsrc/materials/ViscoelasticMaxwell.hh
 *
 * @brief Class for basic Maxwell viscoelastic functions.
 */

#if !defined(pylith_materials_viscoelasticmaxwell_hh)
#define pylith_materials_viscoelasticmaxwell_hh

// Include directives ---------------------------------------------------
#include "materialsfwd.hh" // forward declarations

// ViscoelasticMaxwell --------------------------------------------------
/** @brief Class for basic Maxwell viscoelastic functions.
 *
 * This class contains functions that can be used by any linear
 * Maxwell viscoelastic class.
 */
class pylith::materials::ViscoelasticMaxwell
{ // class ViscoelasticMaxwell

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Compute viscous strain parameter.
   *
   * @param dt Time step.
   * @param maxwellTime Maxwell time.
   *
   * @returns Viscous strain parameter.
   */
  static double viscousStrainParam(const double dt,
				   const double maxwellTime);

}; // class ViscoelasticMaxwell

#endif // pylith_materials_viscoelasticmaxwell_hh


// End of file 
