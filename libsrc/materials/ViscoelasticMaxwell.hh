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
 * @brief C++ ViscoelasticMaxwell object.
 *
 * This class contains a single function that can be used by any
 * linear Maxwell viscoelastic class.
 */

#if !defined(pylith_materials_viscoelasticmaxwell_hh)
#define pylith_materials_viscoelasticmaxwell_hh

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class ViscoelasticMaxwell;
  } // materials

} // pylith

/// C++ abstract base class for Material object.
class pylith::materials::ViscoelasticMaxwell
{ // class ViscoelasticMaxwell

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Compute viscous strain parameter.
   *
   * @returns Viscous strain parameter.
   */
  static double computeVisStrain(const double dt,
				 const double maxwelltime);

}; // class ViscoelasticMaxwell

#endif // pylith_materials_viscoelasticmaxwell_hh


// End of file 
