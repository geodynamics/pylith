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
#include "pylith/materials/AuxiliaryFactoryElasticity.hh" // ISA AuxiliaryFactoryElasticity

class pylith::materials::AuxiliaryFactoryElastic : public pylith::materials::AuxiliaryFactoryElasticity {
    friend class TestAuxiliaryFactoryElastic; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryElastic(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryElastic(void);

    /// Add shear modulus subfield to auxiliary subfields.
    void addShearModulus(void);

    /// Add bulk subfield to auxiliary subfields.
    void addBulkModulus(void);

    /// Add reference stress subfield to auxiliary subfields.
    void addReferenceStress(void);

    /// Add reference strain subfield to auxiliary subfields.
    void addReferenceStrain(void);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryElastic(const AuxiliaryFactoryElastic &); ///< Not implemented.
    const AuxiliaryFactoryElastic& operator=(const AuxiliaryFactoryElastic&); ///< Not implemented

}; // class AuxiliaryFactoryElastic

// End of file
