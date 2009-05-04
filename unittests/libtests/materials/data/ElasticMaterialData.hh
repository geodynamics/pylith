
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

  double dtStableImplicit; ///< Stable time step for implicit time stepping.
  double* density; ///< Density at location.
  double* strain; ///< Strain at location.
  double* stress; ///< Stress at location.
  double* elasticConsts; ///< Elastic constants at location.

  double* initialStress; ///< Initial stress at location.
  double* initialStrain; ///< Initial strain at location.

  double* stateVarsUpdated; ///< Updated state variables at location.

};

#endif // pylith_materials_elasticmaterialdata_hh


// End of file
