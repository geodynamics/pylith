
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

#if !defined(pylith_materials_elasticisotropic3ddata_hh)
#define pylith_materials_elasticisotropic3ddata_hh

#include "MaterialData.hh"

namespace pylith {
  namespace materials {
     class ElasticIsotropic3DData;
  } // pylith
} // materials

class pylith::materials::ElasticIsotropic3DData : public MaterialData
{

public: 

  /// Constructor
  ElasticIsotropic3DData(void);

  /// Destructor
  ~ElasticIsotropic3DData(void);

private:

  static const int _numDBValues;
  static const int _numParameters;
  static const char* _dbValues[];
  static const char* _parameterNames[];

};

#endif // pylith_materials_elasticisotropic3ddata_hh


// End of file
