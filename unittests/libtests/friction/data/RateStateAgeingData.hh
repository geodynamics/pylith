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

  static const double _lengthScale;

  static const double _timeScale;

  static const double _pressureScale;

  static const double _densityScale;

  static const double _dt;

  static const int _numPropertyValues[];

  static const int _numStateVarValues[];

  static const char* _dbPropertyValues[];

  static const char* _dbStateVarValues[];

  static const double _dbProperties[];

  static const double _dbStateVars[];

  static const double _properties[];

  static const double _stateVars[];

  static const double _propertiesNondim[];

  static const double _stateVarsNondim[];

  static const double _friction[];

  static const double _slip[];

  static const double _slipRate[];

  static const double _normalTraction[];

  static const double _stateVarsUpdated[];

};

#endif // pylith_friction_ratestateageingdata_hh

// End of file
