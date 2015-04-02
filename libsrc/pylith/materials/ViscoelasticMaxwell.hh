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

/** @file libsrc/materials/ViscoelasticMaxwell.hh
 *
 * @brief Class for basic Maxwell viscoelastic functions.
 */

#if !defined(pylith_materials_viscoelasticmaxwell_hh)
#define pylith_materials_viscoelasticmaxwell_hh

// Include directives ---------------------------------------------------
#include "materialsfwd.hh" // forward declarations

#include "pylith/utils/types.hh" // USES PylithScalar

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
  static PylithScalar viscousStrainParam(const PylithScalar dt,
				   const PylithScalar maxwellTime);

}; // class ViscoelasticMaxwell

#endif // pylith_materials_viscoelasticmaxwell_hh


// End of file 
