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

  int numDBValues; ///< Number of database values
  int numParameters; ///< Number of parameters

  char** dbValues; ///< Aray of names of database values;
  char** parameterNames; //< Array of names of parameters

  double* dbData; ///< Array of database values
  double* parameterData; ///< Array of parameter values
};

#endif // pylith_materials_materialdata_hh

// End of file
