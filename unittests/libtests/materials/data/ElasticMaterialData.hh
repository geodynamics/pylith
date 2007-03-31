
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

#if !defined(pylith_materials_elasticmaterialdata_hh)
#define pylith_materials_elasticmaterialdata_hh

#include "MaterialData.hh"

namespace pylith {
  namespace materials {
     class ElasticMaterialData;
  } // pylith
} // materials

class pylith::materials::ElasticMaterialData : public MaterialData
{

public: 

  /// Constructor
  ElasticMaterialData(void);

  /// Destructor
  ~ElasticMaterialData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  int numLocs; ///< Number of locations

  double* density; ///< Density at locations
  double* strain; ///< Strain at locations
  double* stress; ///< Stress at locations
  double* elasticConsts; ///< Elastic constants at locations

};

#endif // pylith_materials_elasticmaterialdata_hh


// End of file
