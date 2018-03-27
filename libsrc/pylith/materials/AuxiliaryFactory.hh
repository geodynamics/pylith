// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/AuxiliaryFactory.hh
 *
 * @brief C++ helper class for setting up auxiliary fields for materials.
 */

#if !defined(pylith_materials_auxiliaryfactory_hh)
#define pylith_materials_auxiliaryfactory_hh

// Include directives ---------------------------------------------------
#include "materialsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

// AuxiliaryFactory-----------------------------------------------
/// @brief C++ helper class for setting up auxiliary fields for materials.
class pylith::materials::AuxiliaryFactory : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactory;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactory(void);

    /// Destructor.
    ~AuxiliaryFactory(void);

    /// Add density subfield to auxiliary fields.
    void density(void);

    /// Add shear modulus subfield to auxiliary fields.
    void shearModulus(void);

    /// Add bulk subfield to auxiliary fields.
    void bulkModulus(void);

    /** Add gravity subfield to auxiliary fields.
     *
     * @param[in] gf Gravity field.
     */
    void gravityField(spatialdata::spatialdb::GravityField* gf);

    /// Add body force subfield to auxiliary fields.
    void bodyForce(void);

    /// Add reference stress subfield to auxiliary fields.
    void referenceStress(void);

    /// Add reference strain subfield to auxiliary fields.
    void referenceStrain(void);

  /** Add Maxwell time subfield to auxiliary fields.
   *
   * The identifier is used to distinguish among multiple subfields of
   * the same parameters for cases like a generalized Maxwell model.
   *
   * @param[in] identifier Identifier appended to field name.
   */
    void maxwellTime(const char* identifier =NULL);

    /// Add total strain subfield to auxiliary fields.
    void totalStrain(void);

  /** Add viscous strain subfield to auxiliary fields.
   *
   * The identifier is used to distinguish among multiple subfields of
   * the same parameters for cases like a generalized Maxwell model.
   *
   * @param[in] identifier Identifier appended to field name.
   */
    void viscousStrain(const char* identigier =NULL);

  /** Add shear modulus ratio subfield to auxliary fields.
   *
   * The identifier is used to distinguish among multiple subfields of
   * the same parameters for cases like a generalized Maxwell model.
   *
   * @param[in] identifier Identifier appended to field name.
   */
  void shearModulusRatio(const char* identifier =NULL);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    static const char* _genericComponent; ///< Name of generic component.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    AuxiliaryFactory(const AuxiliaryFactory &);   ///< Not implemented.
    const AuxiliaryFactory& operator=(const AuxiliaryFactory&);   ///< Not implemented

}; // class AuxiliaryFactory

#endif // pylith_materials_auxiliaryfactory_hh


// End of file
