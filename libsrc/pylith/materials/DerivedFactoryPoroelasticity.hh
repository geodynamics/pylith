// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/materials/materialsfwd.hh" // forward declarations
#include "pylith/materials/DerivedFactoryElasticity.hh" // ISA DerivedFactoryElasticity

class pylith::materials::DerivedFactoryPoroelasticity : public pylith::materials::DerivedFactoryElasticity {
    friend class TestDerivedFactoryElasticity; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    DerivedFactoryPoroelasticity(void);

    /// Destructor.
    virtual ~DerivedFactoryPoroelasticity(void);

    /// Add bulk density subfield to derived field
    void addBulkDensity(void);

    /// Add water content subfield to derived field
    void addWaterContent(void);

    /// Add subfields using discretizations provided.
    void addSubfields(void);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    DerivedFactoryPoroelasticity(const DerivedFactoryPoroelasticity &); ///< Not implemented.
    const DerivedFactoryPoroelasticity& operator=(const DerivedFactoryPoroelasticity&); ///< Not implemented

}; // class DerivedFactoryPoroelasticity

// End of file
