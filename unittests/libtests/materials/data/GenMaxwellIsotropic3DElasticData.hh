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

// DO NOT EDIT THIS FILE
// This file was generated from python application genmaxwellisotropic3delastic.

#if !defined(pylith_materials_genmaxwellisotropic3delasticdata_hh)
#define pylith_materials_genmaxwellisotropic3delasticdata_hh

#include "ElasticMaterialData.hh"

namespace pylith {
  namespace materials {
     class GenMaxwellIsotropic3DElasticData;
  } // pylith
} // materials

class pylith::materials::GenMaxwellIsotropic3DElasticData : public ElasticMaterialData
{

public: 

  /// Constructor
  GenMaxwellIsotropic3DElasticData(void);

  /// Destructor
  ~GenMaxwellIsotropic3DElasticData(void);

private:

  static const int _dimension;

  static const int _numDBValues;

  static const int _numParameters;

  static const int _numParamsQuadPt;

  static const int _numLocs;

  static const int _numParamValues[];

  static const char* _dbValues[];

  static const double _dbData[];

  static const double _parameterData[];

  static const double _density[];

  static const double _strain[];

  static const double _stress[];

  static const double _elasticConsts[];

};

#endif // pylith_materials_genmaxwellisotropic3delasticdata_hh

// End of file
