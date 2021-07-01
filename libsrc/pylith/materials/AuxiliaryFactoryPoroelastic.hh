// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/AuxiliaryFactoryElastic.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for elastic materials.
 */

#if !defined(pylith_materials_auxiliaryfactoryporoelastic_hh)
#define pylith_materials_auxiliaryfactoryporoelastic_hh

#include "materialsfwd.hh" // forward declarations
#include "pylith/materials/AuxiliaryFactoryPoroelasticity.hh" // ISA AuxiliaryFactoryPoroelasticity

class pylith::materials::AuxiliaryFactoryPoroelastic : public pylith::materials::AuxiliaryFactoryPoroelasticity {
    friend class TestAuxiliaryFactoryPoroelastic; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryPoroelastic(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryPoroelastic(void);

    /// Add isotropic permeability subfield to auxiliary subfields.
    void addIsotropicPermeability(void);

    /// Add tensor permeability subfield to auxiliary subfields.
    void addTensorPermeability(void);

    /// Add drained Bulk Modulus subfield to auxiliary subfields.
    void addDrainedBulkModulus(void);

    /// Add undrained Bulk Modulus subfield to auxiliary subfields.
    void addUndrainedBulkModulus(void);

    /// Add fluid Bulk Modulus subfield to auxiliary subfields.
    void addFluidBulkModulus(void);

    /// Add fluid Biot Coefficient subfield to auxiliary subfields.
    void addBiotCoefficient(void);

    /// Add fluid Biot Modulus subfield to auxiliary subfields.
    void addBiotModulus(void);

    /// Add reference stress subfield to auxiliary fields.
    void addReferenceStress(void);

    /// Add reference strain subfield to auxiliary fields.
    void addReferenceStrain(void);

    /// Add shear modulus subfield to auxiliary subfields.
    void addShearModulus(void);

    /// Add solid bulk modulus subfield to auxiliary subfields.
    void addSolidBulkModulus(void);

    /// Add shear modulus subfield to auxiliary subfields.
    void addYoungsModulus(void);

    /// Add solid bulk modulus subfield to auxiliary subfields.
    void addPoissonsRatio(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryPoroelastic(const AuxiliaryFactoryPoroelastic &); ///< Not implemented.
    const AuxiliaryFactoryPoroelastic& operator=(const AuxiliaryFactoryPoroelastic&); ///< Not implemented

}; // class AuxiliaryFactoryPoroelastic

#endif // pylith_materials_auxiliaryfactoryporoelastic_hh

// End of file
