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

#if !defined(pylith_materials_materialdata_hh)
#define pylith_materials_materialdata_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
  namespace materials {
     class MaterialData;
  } // pylith
} // materials

class pylith::materials::MaterialData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  MaterialData(void);

  /// Destructor
  ~MaterialData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  int dimension; ///< Number of dimensions
  int numLocs; ///< Number of locations

  int numProperties; ///< Number of parameters for physical properties.
  int numStateVars; ///< Number of state variables.
  int numDBProperties; ///< Number of db values for properties.
  int numDBStateVars; ///< Number of db values for state variables.

  int numPropsQuadPt; ///< Number of properties at each quadrature point
  int numVarsQuadPt; ///< Number of state variables at each quadrature point

  int* numPropertyValues; ///< Number of values for each property
  int* numStateVarValues; ///< Number of values for each state variable

  char** dbPropertyValues; ///< Names of db values for properties.
  char** dbStateVarValues; ///< Names of db values for state variables.

  PylithScalar* dbProperties; ///< Database values for properties at locations.
  PylithScalar* dbStateVars; ///< Database values for state variables at locations.
  PylithScalar* properties; ///< Properties at locations.
  PylithScalar* stateVars; ///< State variables at locations.
  PylithScalar* propertiesNondim; ///< Nondimensional properties at locations.
  PylithScalar* stateVarsNondim; ///< Nondimensional state variables at locations.

  PylithScalar lengthScale; ///< Length scale for nondimensionalization.
  PylithScalar timeScale; ///< Time scale for nondimensionalization.
  PylithScalar pressureScale; ///< Pressure scale for nondimensionalization.
  PylithScalar densityScale; ///< Density scale for nondimensionalization.

};

#endif // pylith_materials_materialdata_hh

// End of file
