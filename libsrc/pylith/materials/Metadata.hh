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

#include <string> // HASA std::string

// MaterialMetadata -----------------------------------------------------
/// @brief C++ object for material metadata.
class pylith::materials::Metadata
{ // Mesh
  friend class TestMetadata; // unit testing

// PUBLIC STRUCTS ///////////////////////////////////////////////////////
public :

  struct ParamDescription {
    std::string name;
    int fiberDim;
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
  
  /** Get number of physical properties.
   *
   * @returns Number of physical properties.
   */
  int numProperties(void) const;

  /** Get physical property index.
   *
   * @param index Physical property index.
   */
  const ParamDescription& getProperty(const int index) const;

  /** Get number of state variables.
   *
   * @returns Number of state variables.
   */
  int numStateVars(void) const;

  /** Get state variable index.
   *
   * @param index State variable index.
   */
  const ParamDescription& getStateVar(const int index) const;

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


// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  ParamDescription* _properties; ///< Physical properties information.
  ParamDescription* _stateVars; ///< State variable information.

  const char* const* _dbProperties; ///< Names of db values for properties. 
  const char* const* _dbStateVars; ///< Names of db values for state varaibles.

  const int _numProperties; ///< Number of physical properties.
  const int _numStateVars; ///< NUmber of state variables.

  const int _numDBProperties; ///< Number of db values for properties.
  const int _numDBStateVars; ///< Number of db values for state variables.
  
// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  const Metadata& operator=(const Metadata&); ///< Not implemented

}; // Metadata

#include "Metadata.icc"

#endif // pylith_materials_metadata_hh


// End of file
