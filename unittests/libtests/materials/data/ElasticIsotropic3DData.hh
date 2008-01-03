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
// This file was generated from python application elasticisotropic3d.

#if !defined(pylith_materials_elasticisotropic3ddata_hh)
#define pylith_materials_elasticisotropic3ddata_hh

#include "ElasticMaterialData.hh"

namespace pylith {
  namespace materials {
     class ElasticIsotropic3DData;
  } // pylith
} // materials

class pylith::materials::ElasticIsotropic3DData : public ElasticMaterialData
{

public: 

  /// Constructor
  ElasticIsotropic3DData(void);

  /// Destructor
  ~ElasticIsotropic3DData(void);

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

#endif // pylith_materials_elasticisotropic3ddata_hh

// End of file
