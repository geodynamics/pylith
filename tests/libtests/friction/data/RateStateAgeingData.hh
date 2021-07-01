// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#if !defined(pylith_friction_ratestateageingdata_hh)
#define pylith_friction_ratestateageingdata_hh

#include "FrictionModelData.hh"

namespace pylith {
  namespace friction {
     class RateStateAgeingData;
  } // pylith
} // friction

class pylith::friction::RateStateAgeingData : public FrictionModelData
{

public: 

  /// Constructor
  RateStateAgeingData(void);

  /// Destructor
  ~RateStateAgeingData(void);

private:

  static const int _numLocs;

  static const int _numProperties;

  static const int _numStateVars;

  static const int _numDBProperties;

  static const int _numDBStateVars;

  static const int _numPropsVertex;

  static const int _numVarsVertex;

  static const PylithScalar _lengthScale;

  static const PylithScalar _timeScale;

  static const PylithScalar _pressureScale;

  static const PylithScalar _densityScale;

  static const PylithScalar _dt;

  static const int _numPropertyValues[];

  static const int _numStateVarValues[];

  static const char* _dbPropertyValues[];

  static const char* _dbStateVarValues[];

  static const PylithScalar _dbProperties[];

  static const PylithScalar _dbStateVars[];

  static const PylithScalar _properties[];

  static const PylithScalar _stateVars[];

  static const PylithScalar _propertiesNondim[];

  static const PylithScalar _stateVarsNondim[];

  static const PylithScalar _friction[];

  static const PylithScalar _frictionDeriv[];

  static const PylithScalar _slip[];

  static const PylithScalar _slipRate[];

  static const PylithScalar _normalTraction[];

  static const PylithScalar _stateVarsUpdated[];

};

#endif // pylith_friction_ratestateageingdata_hh

// End of file
