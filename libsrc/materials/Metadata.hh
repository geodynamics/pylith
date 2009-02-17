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

/**
 * @file libsrc/topology/Metadata.hh
 *
 * @brief C++ object for material metadata.
 *
 * Extends Sieve mesh to include coordinate system associated with
 * domain.
 */

#if !defined(pylith_materials_metadata_hh)
#define pylith_materials_metadata_hh

// Include directives ---------------------------------------------------
#include "materialsfwd.hh" // forward declarations

#include "pylith/topology/FieldBase.hh" // USES FieldBase::VectorFieldEnum
#include "pylith/utils/array.hh" // HASA string_vector

#include <map> // HASA std::map

// MaterialMetadata -----------------------------------------------------
class pylith::materials::Metadata
{ // Mesh
  friend class TestMetadata; // unit testing

// PUBLIC ENUMS /////////////////////////////////////////////////////////
public :

  enum ValueEnum {
    PROPERTY=0, ///< Property value.
    STATEVAR=1, ///< State variable value.
  }; // ValueEnum

// PUBLIC STRUCTS ///////////////////////////////////////////////////////
public :

  struct ParamDescription {
    const char* name;
    const int fiberDim;
    topology::FieldBase::VectorFieldEnum fieldType;
  }; // ParamDescription

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /** Constructor.
   *
   * @param numProperties Number of physical properties for material.
   * @param properties Array of property descriptions.
   * @param numDBProperties Number of database values for physical properties.
   * @param dbProperties Array of names for database values for properties.
   * @param numStateVars Number of state variables for material.
   * @param stateVars Array of state variable descriptions.
   * @param numDBStateVars Number of database values for state variables.
   * @param dbStateVars Array of names for database values for state variables.
   */
  Metadata(const int numProperties,
	   const ParamDescription* properties,
	   const int numDBProperties,
	   const char* dbProperties[],
	   const int numStateVars,
	   const ParamDescription* stateVars,
	   const int numDBStateVars,
	   const char* dbStateVars[]);

  /// Default destructor
  ~Metadata(void);

  /** Get names of properties.
   * 
   * @returns Array of names of properties.
   */
  const string_vector& properties(void) const;

  /** Get names of state variables.
   *
   * @returns Array of names of state variables.
   */
  const string_vector& stateVars(void) const;

  /** Get fiber dimension of value.
   *
   * @param name Name of value.
   * @param valueType Type of value.
   */
  int fiberDim(const char* name,
	       const ValueEnum valueType) const;

  /** Get type of vector field associated with value.
   *
   * @param name Name of value.
   * @param valueType Type of value.
   */
  topology::FieldBase::VectorFieldEnum
  fieldType(const char* name,
	    const ValueEnum valueType) const;
   
  /** Get names of database values for physical properties.
   *
   * @returns Array of names.
   */
  const char* const* dbProperties(void) const;

  /** Get number of database values for physical properties.
   *
   * @returns Number of database values.
   */
  int numDBProperties(void) const;

  /** Get names of database values for state variables.
   *
   * @returns Array of names.
   */
  const char* const* dbStateVars(void) const;

  /** Get number of database values for state variables.
   *
   * @returns Number of database values.
   */
  int numDBStateVars(void) const;


// PRIVATE STRUCTS //////////////////////////////////////////////////////
private :

  struct ParameterInfo {
    int fiberDim; ///< Fiber dimension for parameter.
    topology::FieldBase::VectorFieldEnum fieldType; ///< Type of Vector field.
  }; // ParameterInfo

// PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private :

  typedef std::map< std::string, ParameterInfo > ParameterMap;

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  ParameterMap _properties; ///< Physical properties information.
  ParameterMap _stateVars; ///< State variable information.

  string_vector _propertyNames; ///< Names of physical properties.
  string_vector _stateVarNames; ///< Names of state variables.
  
  const char* const* _dbProperties; ///< Names of db values for properties. 
  const char* const* _dbStateVars; ///< Names of db values for state varaibles.
  const int _numDBProperties; ///< Number of db values for properties.
  const int _numDBStateVars; ///< Number of db values for state variables.
  
// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Metadata(const Metadata&); ///< Not implemented
  const Metadata& operator=(const Metadata&); ///< Not implemented

}; // Metadata

#include "Metadata.icc"

#endif // pylith_materials_metadata_hh


// End of file
