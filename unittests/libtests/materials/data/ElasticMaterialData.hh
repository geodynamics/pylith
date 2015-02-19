
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

  PylithScalar dtStableImplicit; ///< Stable time step for implicit time stepping.
  PylithScalar dtStableExplicit; ///< Stable time step for explicit time stepping.
  PylithScalar* density; ///< Density at location.
  PylithScalar* strain; ///< Strain at location.
  PylithScalar* stress; ///< Stress at location.
  PylithScalar* elasticConsts; ///< Elastic constants at location.

  PylithScalar* initialStress; ///< Initial stress at location.
  PylithScalar* initialStrain; ///< Initial strain at location.

  PylithScalar* stateVarsUpdated; ///< Updated state variables at location.

};

#endif // pylith_materials_elasticmaterialdata_hh


// End of file
