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

/** @file modulesrc/materials/ElasticStrain1D.i
 *
 * Python interface to C++ ElasticStrain1D object.
 */

namespace pylith {
  namespace materials {

    class ElasticStrain1D : public ElasticMaterial
    { // class ElasticStrain1D

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor
      ElasticStrain1D(void);

      /// Destructor
      ~ElasticStrain1D(void);

    }; // class ElasticStrain1D

  } // materials
} // pylith


// End of file
