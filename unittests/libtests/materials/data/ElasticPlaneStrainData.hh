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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application elasticplanestrain.

#if !defined(pylith_materials_elasticplanestraindata_hh)
#define pylith_materials_elasticplanestraindata_hh

#include "ElasticMaterialData.hh"

namespace pylith {
  namespace materials {
     class ElasticPlaneStrainData;
  } // pylith
} // materials

class pylith::materials::ElasticPlaneStrainData : public ElasticMaterialData
{

public: 

  /// Constructor
  ElasticPlaneStrainData(void);

  /// Destructor
  ~ElasticPlaneStrainData(void);

private:

  static const int _dimension;

  static const int _numLocs;

  static const int _numProperties;

  static const int _numStateVars;

  static const int _numDBProperties;

  static const int _numDBStateVars;

  static const int _numPropsQuadPt;

  static const int _numVarsQuadPt;

  static const double _lengthScale;

  static const double _timeScale;

  static const double _pressureScale;

  static const double _densityScale;

  static const double _dtStableImplicit;

  static const int _numPropertyValues[];

  static const int* _numStateVarValues;

  static const char* _dbPropertyValues[];

  static const char** _dbStateVarValues;

  static const double _dbProperties[];

  static const double* _dbStateVars;

  static const double _properties[];

  static const double* _stateVars;

  static const double _propertiesNondim[];

  static const double* _stateVarsNondim;

  static const double _density[];

  static const double _strain[];

  static const double _stress[];

  static const double _elasticConsts[];

  static const double _initialStress[];

  static const double _initialStrain[];

  static const double* _stateVarsUpdated;

};

#endif // pylith_materials_elasticplanestraindata_hh

// End of file
