// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/materials/materialsfwd.hh" // forward declarations
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

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryPoroelastic(const AuxiliaryFactoryPoroelastic &); ///< Not implemented.
    const AuxiliaryFactoryPoroelastic& operator=(const AuxiliaryFactoryPoroelastic&); ///< Not implemented

}; // class AuxiliaryFactoryPoroelastic

// End of file
