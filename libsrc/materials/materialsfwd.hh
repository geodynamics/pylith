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
    class GenMaxwellIsotropic3D;

  } // materials
} // pylith


#endif // pylith_materials_materialsfwd_hh


// End of file 
