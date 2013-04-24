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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#if !defined(pylith_feassemble_elasticityexplicitdata_hh)
#define pylith_feassemble_elasticityexplicitdata_hh

#include "IntegratorData.hh" // ISA IntegratorData

namespace pylith {
  namespace feassemble {
     class ElasticityExplicitData;
  } // pylith
} // feassemble

class pylith::feassemble::ElasticityExplicitData : public IntegratorData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  ElasticityExplicitData(void);

  /// Destructor
  ~ElasticityExplicitData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  /// @name Scales information for nondimensionalization.
  //@{
  PylithScalar lengthScale; ///< Length scale.
  PylithScalar pressureScale; ///< Pressure scale.
  PylithScalar timeScale; ///< Time scale.
  PylithScalar densityScale; ///< Density scale.
  //@}  

  /// @name Calculated values.
  //@{
  PylithScalar dtStableExplicit; ///< Stable time step for explicit time integration.
  //@}
};

#endif // pylith_feassemble_elasticityexplicitdata_hh

// End of file
