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
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

// AuxiliaryFactory-----------------------------------------------
/// @brief C++ helper class for setting up auxiliary fields for materials.
class pylith::materials::AuxiliaryFactory : public pylith::utils::GenericComponent {
    friend class TestAuxiliaryFactory;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param[inout material Material for auxiliary fields.
     * @param[in] normalizer Nondimensionalizer for problem.
     * @param[in] spaceDim Spatial dimension of problem.
     */
    AuxiliaryFactory(const MaterialNew& material,
                     const spatialdata::units::Nondimensional& normalizer,
                     const int spaceDim);

    /// Destructor.
    ~AuxiliaryFactory(void);

    /// Add density field to auxiliary fields.
    void density(void) const;

    /// Add shear modulus field to auxiliary fields.
    void shearModulus(void) const;

    /// Add bulk field to auxiliary fields.
    void bulkModulus(void) const;

    /** Add gravity field to auxiliary fields.
     *
     * @param[in] gf Gravity field.
     */
    void gravityField(spatialdata::spatialdb::GravityField* gf) const;

    /// Add body force field to auxiliary fields.
    void bodyForce(void) const;

    /// Add reference stress field to auxiliary fields.
    void referenceStress(void) const;

    /// Add reference strain field to auxiliary fields.
    void referenceStrain(void) const;

    /// Add Maxwell time field to auxiliary fields.
    void maxwellTime(void) const;

    /// Add total strain field to auxiliary fields.
    void totalStrain(void) const;

    /// Add viscous strain field to auxiliary fields.
    void viscousStrain(void) const;


    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    const MaterialNew& _material; ///< Material with auxiliary fields.
    const spatialdata::units::Nondimensional& _normalizer; ///< Nondimensionalizer.
    const int _spaceDim; ///< Spatal dimension of problem.

    static const char* _genericComponent; ///< Name of generic component.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    AuxiliaryFactory(const AuxiliaryFactory &);   ///< Not implemented.
    const AuxiliaryFactory& operator=(const AuxiliaryFactory&);   ///< Not implemented

}; // class AuxiliaryFactory

#endif // pylith_materials_auxiliaryfactory_hh


// End of file
