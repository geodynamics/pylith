// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/** @file libsrc/materials/materialsfwd.hh
 *
 * @brief Forward declarations for PyLith materials objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_materials_materialsfwd_hh)
#define pylith_materials_materialsfwd_hh

namespace pylith {
  namespace materials {

    class Metadata;
    class ElasticMetadata;

    class Material;

    class ElasticMaterial;
    class ElasticStrain1D;
    class ElasticStress1D;
    class ElasticPlaneStrain;
    class ElasticPlaneStress;
    class ElasticIsotropic3D;
    class MaxwellIsotropic3D;
    class MaxwellPlaneStrain;
    class GenMaxwellIsotropic3D;
    class GenMaxwellPlaneStrain;
    class GenMaxwellQpQsIsotropic3D;
    class PowerLaw3D;
    class PowerLawPlaneStrain;
    class DruckerPrager3D;
    class DruckerPragerPlaneStrain;

    class EffectiveStress;
    class ViscoelasticMaxwell;

  } // materials
} // pylith


#endif // pylith_materials_materialsfwd_hh


// End of file 
