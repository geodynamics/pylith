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

#if !defined(pylith_materials_materialdata_hh)
#define pylith_materials_materialdata_hh

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

  double* dbProperties; ///< Database values for properties at locations.
  double* dbStateVars; ///< Database values for state variables at locations.
  double* properties; ///< Properties at locations.
  double* stateVars; ///< State variables at locations.
  double* propertiesNondim; ///< Nondimensional properties at locations.
  double* stateVarsNondim; ///< Nondimensional state variables at locations.

};

#endif // pylith_materials_materialdata_hh

// End of file
