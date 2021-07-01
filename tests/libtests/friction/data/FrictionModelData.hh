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

#if !defined(pylith_friction_frictionmodeldata_hh)
#define pylith_friction_frictionmodeldata_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
  namespace friction {
     class FrictionModelData;
  } // pylith
} // friction

class pylith::friction::FrictionModelData {

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  FrictionModelData(void);

  /// Destructor
  ~FrictionModelData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  int numLocs; ///< Number of locations

  int numProperties; ///< Number of parameters for physical properties.
  int numStateVars; ///< Number of state variables.
  int numDBProperties; ///< Number of db values for properties.
  int numDBStateVars; ///< Number of db values for state variables.

  int numPropsVertex; ///< Number of properties at each vertex.
  int numVarsVertex; ///< Number of state variables at each vertex.

  int* numPropertyValues; ///< Number of values for each property
  int* numStateVarValues; ///< Number of values for each state variable

  char** dbPropertyValues; ///< Names of db values for properties.
  char** dbStateVarValues; ///< Names of db values for state variables.

  PylithScalar dt; ///< Time step
  PylithScalar* dbProperties; ///< Database values for properties at locations.
  PylithScalar* dbStateVars; ///< Database values for state variables at locations.
  PylithScalar* properties; ///< Properties at locations.
  PylithScalar* stateVars; ///< State variables at locations.
  PylithScalar* propertiesNondim; ///< Nondimensional properties at locations.
  PylithScalar* stateVarsNondim; ///< Nondimensional state variables at locations.

  PylithScalar* friction; ///< Friction at locations.
  PylithScalar* frictionDeriv; ///< Derivative of friction with slip at locations.
  PylithScalar* slip; ///< Slip at locations.
  PylithScalar* slipRate; ///< Slip rate at locations.
  PylithScalar* normalTraction; ///< Normal traction at locations.

  PylithScalar* stateVarsUpdated; ///< Updated state variables at location.

  PylithScalar lengthScale; ///< Length scale for nondimensionalization.
  PylithScalar timeScale; ///< Time scale for nondimensionalization.
  PylithScalar pressureScale; ///< Pressure scale for nondimensionalization.
  PylithScalar densityScale; ///< Density scale for nondimensionalization.

};

#endif // pylith_friction_frictionmodeldata_hh

// End of file
