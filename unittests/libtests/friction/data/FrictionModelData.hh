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

#if !defined(pylith_friction_frictionmodeldata_hh)
#define pylith_friction_frictionmodeldata_hh

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

  double dt; ///< Time step
  double* dbProperties; ///< Database values for properties at locations.
  double* dbStateVars; ///< Database values for state variables at locations.
  double* properties; ///< Properties at locations.
  double* stateVars; ///< State variables at locations.
  double* propertiesNondim; ///< Nondimensional properties at locations.
  double* stateVarsNondim; ///< Nondimensional state variables at locations.

  double* friction; ///< Friction at locations.
  double* slip; ///< Slip at locations.
  double* slipRate; ///< Slip rate at locations.
  double* normalTraction; ///< Normal traction at locations.

  double* stateVarsUpdated; ///< Updated state variables at location.

  double lengthScale; ///< Length scale for nondimensionalization.
  double timeScale; ///< Time scale for nondimensionalization.
  double pressureScale; ///< Pressure scale for nondimensionalization.
  double densityScale; ///< Density scale for nondimensionalization.

};

#endif // pylith_friction_frictionmodeldata_hh

// End of file
