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

  int numDBValues; ///< Number of database values
  int numParameters; ///< Number of parameters
  int numParamsQuadPt; ///< Number of parameters per quadrature point.
  int* numParamValues; ///< Number of values for each parameter

  char** dbValues; ///< Aray of names of database values;

  double* dbData; ///< Array of database values at locations
  double* parameterData; ///< Array of parameter values at locations
};

#endif // pylith_materials_materialdata_hh

// End of file
