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
#include "pylith/materials/AuxiliaryFactoryElastic.hh" // ISA AuxiliaryFactoryElastic

/// @brief C++ helper class for setting up auxiliary fields for materials.
class pylith::materials::AuxiliaryFactoryViscoelastic : public pylith::materials::AuxiliaryFactoryElastic {
    friend class TestAuxiliaryFactoryViscoelastic; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryViscoelastic(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryViscoelastic(void);

    /// Add Maxwell time subfield to auxiliary subfields.
    void addMaxwellTime(void);

    /// Add Maxwell time subfield for Generalized Maxwell to auxiliary subfields.
    void addMaxwellTimeGeneralizedMaxwell(void);

    /// Add shear modulus ratio subfield for Generalized Maxwell to auxiliary subfields.
    void addShearModulusRatioGeneralizedMaxwell(void);

    /// Add power-law reference strain rate subfield to auxiliary subfields.
    void addPowerLawReferenceStrainRate(void);

    /// Add power-law reference stress subfield to auxiliary subfields.
    void addPowerLawReferenceStress(void);

    /// Add power-law exponent subfield to auxiliary subfields.
    void addPowerLawExponent(void);

    /// Add total strain subfield to auxiliary subfields.
    void addTotalStrain(void);

    /// Add stress subfield to auxiliary subfields.
    void addDeviatoricStress(void);

    /// Add viscous strain subfield to auxiliary subfields.
    void addViscousStrain(void);

    /// Add viscous strain subfield for Generalized Maxwell to auxiliary subfields.
    void addViscousStrainGeneralizedMaxwell(void);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryViscoelastic(const AuxiliaryFactoryViscoelastic &); ///< Not implemented.
    const AuxiliaryFactoryViscoelastic& operator=(const AuxiliaryFactoryViscoelastic&); ///< Not implemented

}; // class AuxiliaryFactoryViscoelastic

// End of file
