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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/materials/Metadata.hh
 *
 * @brief C++ object for material metadata.
 */

#if !defined(pylith_materials_metadata_hh)
#define pylith_materials_metadata_hh

// Include directives ---------------------------------------------------
#include "materialsfwd.hh" // forward declarations

#include "pylith/topology/FieldBase.hh" // USES FieldBase::VectorFieldEnum
#include "pylith/utils/array.hh" // HASA string_vector

#include <map> // HASA std::map

// MaterialMetadata -----------------------------------------------------
/** @brief C++ object for material metadata.
 *
 * Extends Sieve mesh to include coordinate system associated with
 * domain.
 */

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
   * @param properties Array of property descriptions.
   * @param numProperties Number of physical properties for material.
   * @param dbProperties Array of names for database values for properties.
   * @param numDBProperties Number of database values for physical properties.
   * @param stateVars Array of state variable descriptions.
   * @param numStateVars Number of state variables for material.
   * @param dbStateVars Array of names for database values for state variables.
   * @param numDBStateVars Number of database values for state variables.
   */
  Metadata(const ParamDescription* properties,
	   const int numProperties,
	   const char* dbProperties[],
	   const int numDBProperties,
	   const ParamDescription* stateVars,
	   const int numStateVars,
	   const char* dbStateVars[],
	   const int numDBStateVars);

  /** Copy constructor.
   *
   * @parameter m Metadataw to copy.
   */
  Metadata(const Metadata& m);

  /// Default destructor
  ~Metadata(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
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

  const Metadata& operator=(const Metadata&); ///< Not implemented

}; // Metadata

#include "Metadata.icc"

#endif // pylith_materials_metadata_hh


// End of file
